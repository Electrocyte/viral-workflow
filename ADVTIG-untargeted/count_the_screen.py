#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import os, re
import glob, json
import sys
import pandas as pd
from multiprocessing import Pool, cpu_count
from typing import Dict, Tuple, Optional


# Filtering knobs
TP_REQUIRE_PRIMARY = True          # drop tp:A:S
# NM_RATIO_THRESH    = 0.10          # edit-distance threshold
# NM_RATIO_THRESH    = 0.15          # edit-distance threshold
NM_DENOM           = "aln"         # "aln" (recommended) or "query"
MIN_ALN_LEN        = 0             # optional


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Count the counter screen results and generate summary statistics.")
    parser.add_argument("-o", "--out-dir", type=str, required=True, help="Directory to save alignment outputs.")
    parser.add_argument("-s", "--samples-config", required=True, help="Sample configuration CSV (same as glue.py)")
    parser.add_argument("--skip", action="store_true", help="Skip existing outputs if present.")
    parser.add_argument("-e", "--edit-distance-threshold", type=float, default=0.1, help="Edit distance threshold for filtering alignments.")    

    return parser.parse_args()


def _hit_key(mapq: int, nm_ratio: float, aln_len: int) -> Tuple[int, float, int]:
    """
    Ranking for "best hit per read":
      1) higher MAPQ
      2) lower NM ratio
      3) higher aligned length
    """
    return (mapq, -nm_ratio, aln_len)


def _parse_parent_from_filename(paf_path: str) -> Optional[Tuple[str, str]]:
    """
    Extract (SeqID, TaxID) from basename like:
      AF288220.1_taxID-9606_mammalian.paf
    """
    base = os.path.basename(paf_path)
    m = re.match(r"^(?P<seqid>.+)_taxID-(?P<taxid>\d+)_mammalian\.paf$", base)
    if not m:
        return None
    return m.group("seqid"), m.group("taxid")


def _infer_sample_from_path(paf_path: str, counter_root: str) -> str:
    """
    Infer sample as the first directory below COUNTER_SCREEN_OUT:
      COUNTER_SCREEN_OUT/<sample>/mammalian/<file>.paf
    """
    rel = os.path.relpath(paf_path, counter_root)
    parts = rel.split(os.sep)
    return parts[0] if parts else ""


def _parse_paf_tags(tag_fields) -> Dict[str, str]:
    """Extract only tags we care about from PAF optional fields."""
    out = {}
    for t in tag_fields:
        # Examples: "tp:A:P", "NM:i:15", "AS:i:422"
        if t.startswith("tp:A:"):
            out["tp"] = t.split(":")[-1]
        elif t.startswith("NM:i:"):
            out["NM"] = t.split(":")[-1]
        elif t.startswith("AS:i:"):
            out["AS"] = t.split(":")[-1]
    return out


def parse_parent_paf_besthits(
    paf_path: str,
    nm_ratio_thresh: float,
    tp_require_primary: bool = True,
    nm_denom: str = "aln",   # "aln" or "query"
    min_aln_len: int = 0,
) -> Tuple[Dict[str, Dict], int]:
    """
    Stream-parse a parent PAF and select *one best hit per read*.

    Returns:
      best_by_read: read_id -> {"target": accession, "mapq": int, "nm_ratio": float, "aln_len": int}
      seen_reads:  number of unique read_ids that appeared in this PAF (aligned-to-something-on-panel)
    """
    best_by_read: Dict[str, Dict] = {}
    seen = set()

    if (not os.path.exists(paf_path)) or os.stat(paf_path).st_size == 0:
        return best_by_read, 0

    with open(paf_path, "r") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            read_id = fields[0]
            seen.add(read_id)

            try:
                qlen = int(fields[1])
                aln_len = int(fields[10])  # alignment block length
                mapq = int(fields[11])
            except ValueError:
                continue

            if aln_len < min_aln_len:
                continue

            target = fields[5]

            tags = _parse_paf_tags(fields[12:])
            if tp_require_primary:
                tp = tags.get("tp")
                if tp != "P":  # only keep primary
                    continue

            nm_s = tags.get("NM")
            if nm_s is None:
                continue
            try:
                nm = int(nm_s)
            except ValueError:
                continue

            denom = aln_len if nm_denom == "aln" else qlen
            if denom <= 0:
                continue

            nm_ratio = nm / float(denom)
            if nm_ratio > nm_ratio_thresh:
                continue

            new_key = _hit_key(mapq, nm_ratio, aln_len)

            prev = best_by_read.get(read_id)
            if prev is None:
                best_by_read[read_id] = {
                    "target": target,
                    "mapq": mapq,
                    "nm_ratio": nm_ratio,
                    "aln_len": aln_len,
                }
            else:
                prev_key = _hit_key(prev["mapq"], prev["nm_ratio"], prev["aln_len"])
                if new_key > prev_key:
                    best_by_read[read_id] = {
                        "target": target,
                        "mapq": mapq,
                        "nm_ratio": nm_ratio,
                        "aln_len": aln_len,
                    }

    return best_by_read, len(seen)


def load_parent_stats(stats_glob: str) -> pd.DataFrame:

    rows = []
    for p in glob.glob(stats_glob, recursive=True):
        if not os.path.exists(p):
            continue
        try:
            df = pd.read_csv(p)
        except Exception:
            continue
        if df.empty:
            continue
        rows.append(df)
    if not rows:
        return pd.DataFrame()

    out = pd.concat(rows, ignore_index=True)

    # normalize key types for merging
    for c in ["sample", "seqid", "taxid"]:
        if c in out.columns:
            out[c] = out[c].astype(str)

    return out


def build_parent_summary_from_pafs(
    counter_screen_out: str,
    paf_glob: str,
    nm_ratio_thresh: float,
    stats_df: Optional[pd.DataFrame] = None,
    tp_require_primary: bool = True,
    nm_denom: str = "aln",
    min_aln_len: int = 0,
) -> pd.DataFrame:
    """
    Build a per-parent summary table:
      one row per (Sample, SeqID, TaxID) parent PAF file.
    """
    pafs = glob.glob(paf_glob, recursive=True)

    out_rows = []
    for paf in pafs:
        parent = _parse_parent_from_filename(paf)
        if parent is None:
            continue
        seqid, taxid = parent
        sample = _infer_sample_from_path(paf, counter_screen_out)

        best_by_read, seen_reads = parse_parent_paf_besthits(
            paf,
            nm_ratio_thresh=nm_ratio_thresh,
            tp_require_primary=tp_require_primary,
            nm_denom=nm_denom,
            min_aln_len=min_aln_len,
        )

        # count reads per target accession using the BEST hit per read
        acc_to_reads: Dict[str, set] = {}
        for rid, info in best_by_read.items():
            acc = info["target"]
            acc_to_reads.setdefault(acc, set()).add(rid)

        counts = {acc: len(rset) for acc, rset in acc_to_reads.items()}
        passing_any = sum(counts.values())  # unique reads with a kept best hit (should equal len(best_by_read))

        # compute pct if we can
        denom_reads = None
        if stats_df is not None and not stats_df.empty:
            # Try to find a matching stats row (if you generated stats alongside PAFs)
            m = stats_df[
                (stats_df["sample"].astype(str) == str(sample)) &
                (stats_df["seqid"].astype(str) == str(seqid)) &
                (stats_df["taxid"].astype(str) == str(taxid))
            ]
            if not m.empty:
                # If you later add total_reads column, prefer it.
                if "total_reads" in m.columns:
                    try:
                        denom_reads = int(m.iloc[0]["total_reads"])
                    except Exception:
                        denom_reads = None
                else:
                    # Fallback: cannot reliably infer total reads from rounded pct.
                    denom_reads = None

        # If denom unknown, fall back to "seen in PAF" (aligned-to-panel) so pct is at least defined.
        # This is NOT the same as original subset size; it only counts reads that aligned somewhere.
        if denom_reads is None:
            denom_reads = seen_reads

        passing_pct = round((passing_any / float(denom_reads) * 100.0), 2) if denom_reads else 0.0
        has_hits = passing_any > 0

        # top accession
        top_acc = ""
        top_n = 0
        if counts:
            top_acc, top_n = max(counts.items(), key=lambda kv: kv[1])
        top_pct = round((top_n / float(denom_reads) * 100.0), 2) if denom_reads else 0.0

        paf = paf.split("/")[-1].replace(".paf", "")  # relative path

        out_rows.append({
            "Sample": str(sample),
            "SeqID": str(seqid),
            "TaxID": str(taxid),
            "refseq_has_hits": bool(has_hits),
            "refseq_passing_reads": int(passing_any),
            "refseq_passing_pct": float(passing_pct),
            "refseq_top_accession": str(top_acc),
            "refseq_top_reads": int(top_n),
            "refseq_top_pct": float(top_pct),
            "refseq_accession_counts_json": json.dumps(counts, sort_keys=True),
            "refseq_paf_path": paf,
        })

    return pd.DataFrame(out_rows)


def merge_into_df_filter(df_filter: pd.DataFrame, parent_summary: pd.DataFrame) -> pd.DataFrame:
    """
    Merge summary columns into df_filter using (Sample, SeqID, TaxID).
    Handles dtype normalization and fills missing values.
    """
    out = df_filter.copy()

    # normalize types on left
    for c in ["Sample", "SeqID", "TaxID"]:
        if c in out.columns:
            out[c] = out[c].astype(str)

    # normalize on right
    for c in ["Sample", "SeqID", "TaxID"]:
        if c in parent_summary.columns:
            parent_summary[c] = parent_summary[c].astype(str)

    merged = out.merge(
        parent_summary.drop(columns=["refseq_paf_path"], errors="ignore"),
        on=["Sample", "SeqID", "TaxID"],
        how="left",
    )

    # fill defaults
    if "refseq_has_hits" in merged.columns:
        merged["refseq_has_hits"] = merged["refseq_has_hits"].fillna(False)
    for c in ["refseq_passing_reads", "refseq_top_reads"]:
        if c in merged.columns:
            merged[c] = merged[c].fillna(0).astype(int)
    for c in ["refseq_passing_pct", "refseq_top_pct"]:
        if c in merged.columns:
            merged[c] = merged[c].fillna(0.0)
    for c in ["refseq_top_accession", "refseq_accession_counts_json"]:
        if c in merged.columns:
            merged[c] = merged[c].fillna("")

    return merged


def build_sample_rollup(parent_summary: pd.DataFrame) -> pd.DataFrame:
    """
    Expand accession counts JSON and aggregate per sample:
      Sample, accession, reads, n_parents
    """
    if parent_summary.empty:
        return pd.DataFrame(columns=["Sample", "accession", "reads", "n_parents"])

    rows = []
    for _, r in parent_summary.iterrows():
        sample = r["Sample"]
        try:
            d = json.loads(r.get("refseq_accession_counts_json", "{}") or "{}")
        except Exception:
            d = {}
        for acc, n in d.items():
            rows.append({"Sample": sample, "accession": acc, "reads": int(n), "parent": f"{r['SeqID']}_taxID-{r['TaxID']}"})

    if not rows:
        return pd.DataFrame(columns=["Sample", "accession", "reads", "n_parents"])

    df = pd.DataFrame(rows)
    roll = df.groupby(["Sample", "accession"], as_index=False).agg(
        reads=("reads", "sum"),
        n_parents=("parent", pd.Series.nunique),
    )
    return roll.sort_values(["Sample", "reads"], ascending=[True, False])


def build_sample_unique_rollup_from_pafs(
    sample: str,
    paf_glob: str,
    nm_ratio_thresh: float,    
    tp_require_primary: bool = True,
    nm_denom: str = "aln",
    min_aln_len: int = 0,
) -> pd.DataFrame:
    """
    Unique read rollup at sample level:
      one row per (Sample, accession) with unique read counts.
    """
    pafs = glob.glob(paf_glob, recursive=True)

    best_global: Dict[str, Dict] = {}  # read_id -> {"target","mapq","nm_ratio","aln_len"}

    for paf in pafs:
        if not os.path.exists(paf) or os.stat(paf).st_size == 0:
            continue

        with open(paf, "r") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue

                read_id = fields[0]
                try:
                    qlen = int(fields[1])
                    aln_len = int(fields[10])
                    mapq = int(fields[11])
                except ValueError:
                    continue
                if aln_len < min_aln_len:
                    continue

                target = fields[5]
                tags = _parse_paf_tags(fields[12:])

                if tp_require_primary and tags.get("tp") != "P":
                    continue

                nm_s = tags.get("NM")
                if nm_s is None:
                    continue
                try:
                    nm = int(nm_s)
                except ValueError:
                    continue

                denom = aln_len if nm_denom == "aln" else qlen
                if denom <= 0:
                    continue
                nm_ratio = nm / float(denom)
                if nm_ratio > nm_ratio_thresh:
                    continue

                new_key = _hit_key(mapq, nm_ratio, aln_len)

                prev = best_global.get(read_id)
                if prev is None:
                    best_global[read_id] = {"target": target, "mapq": mapq, "nm_ratio": nm_ratio, "aln_len": aln_len}
                else:
                    prev_key = _hit_key(prev["mapq"], prev["nm_ratio"], prev["aln_len"])
                    if new_key > prev_key:
                        best_global[read_id] = {"target": target, "mapq": mapq, "nm_ratio": nm_ratio, "aln_len": aln_len}

    # Unique read counts per accession
    acc_counts: Dict[str, int] = {}
    for rid, info in best_global.items():
        acc = info["target"]
        acc_counts[acc] = acc_counts.get(acc, 0) + 1

    out = pd.DataFrame(
        [{"Sample": sample, "accession": acc, "unique_reads": n} for acc, n in acc_counts.items()]
    ).sort_values(["unique_reads"], ascending=[False])

    return out


def load_samples(sample_config):
    """Load sample names from configuration file."""
    try:
        samples = pd.read_csv(sample_config)

        # Create sample name from components (matching your existing pattern)
        samples["sample"] = samples["date"].astype("str") + "_" + \
            samples["NA"] + "_" + \
            samples["strain"] + "_" + \
            samples["concentration_CFU"] + "_" + \
            samples["batch"].astype("str") + "_" + \
            samples["duration_h"].astype("str")

        sample_names = samples["sample"].unique()

        # Also identify negative control samples
        negative_df = samples.loc[samples["Spike_species"].str.contains("Negative control")]
        negative_samples = list(negative_df["sample"].unique())
        other_samples = [sample for sample in sample_names if sample not in negative_samples]

        return sample_names, negative_samples, other_samples
    except Exception as e:
        print(f"Error loading sample configuration: {e}")
        sys.exit(1)


def main():
    args = parse_arguments()
    samples_config = args.samples_config
    out_dir = args.out_dir
    skip = args.skip

    NM_RATIO_THRESH    = args.edit_distance_threshold          # edit-distance threshold

    samples, _, _ = load_samples(samples_config)
    
    VIRAL_REFSEQ_TSV = "/mnt/usersData/ADVTIG_fa/viral-input-sequences.tsv"
    BACTERIAL_REFSEQ_TSV="/mnt/usersData/ADVTIG_fa/bacterial-input-sequences.tsv"

    cs_out = f"{out_dir}/counter_screen_refseq/"
    csb_out = f"{out_dir}/counter_screen_bacterial_refseq/"
    m_out = f"{out_dir}/analysis/MetaFilt-Third-Pass-untargeted-uviral25-2-complexity.csv"

    print(m_out)
    df_filter = pd.read_csv(m_out)

    if skip:
        # keys
        keys = ["Sample", "SeqID", "TaxID"]

        # columns to bring from df_b
        bcols = [
            "ADV5", "ADV5_reads", "ADV5_pct",
            "DNA_CS", "DNA_CS_reads", "DNA_CS_pct",
            "Cat", "Cat_reads", "Cat_pct",
            "EBV_B95_8", "EBV_B95_8_reads", "EBV_B95_8_pct",
            "FeLV", "FeLV_reads", "FeLV_pct",
            "FeLV_Kawakami_Theilen", "FeLV_Kawakami_Theilen_reads", "FeLV_Kawakami_Theilen_pct",
            "HCoV_NL63", "HCoV_NL63_reads", "HCoV_NL63_pct",
            "Human", "Human_reads", "Human_pct",
            "Monkey", "Monkey_reads", "Monkey_pct",
            "NC_001510_1", "NC_001510_1_reads", "NC_001510_1_pct",
            "NC_001514_1", "NC_001514_1_reads", "NC_001514_1_pct",
            "PCV1", "PCV1_reads", "PCV1_pct",
            "AF038600_1", "AF038600_1_reads", "AF038600_1_pct",
            "Pig", "Pig_reads", "Pig_pct",
            "REOvirus", "REOvirus_reads", "REOvirus_pct",
            "RSV_A", "RSV_A_reads", "RSV_A_pct",
            "Tamarin", "Tamarin_reads", "Tamarin_pct",
            "mammal_total_reads", "mammal_total_pct",
        ]

        refseq_cols = [
            "refseq_has_hits",
            "refseq_passing_reads",
            "refseq_passing_pct",
            "refseq_top_accession",
            "refseq_top_reads",
            "refseq_top_pct",
            "refseq_accession_counts_json",
            "accession_name",
        ]

        # merge two paths
        # counter_screen_out = f"{out_dir}/counter_screen_comb3/"
        save_summary = m_out.replace(".csv", f"-counter-screened-summary-combined3.csv")  
        viral_MERGED_FILTER  = f"{cs_out}/MetaFilt-Third-Pass-untargeted-uviral25-2-complexity-counter-screen-refseq-viral.csv"  
        bacterial_MERGED_FILTER  = f"{csb_out}/MetaFilt-Third-Pass-untargeted-uviral25-2-complexity-counter-screen-refseq-bacterial.csv"
        
        print(viral_MERGED_FILTER)
        print(bacterial_MERGED_FILTER)
        print(save_summary)

        if os.path.exists(viral_MERGED_FILTER) and os.path.exists(bacterial_MERGED_FILTER) and os.path.exists(save_summary):

            df_a = pd.read_csv(viral_MERGED_FILTER)
            df_b = pd.read_csv(bacterial_MERGED_FILTER)
            df_c = pd.read_csv(save_summary)
            print(df_a.columns)
            print(df_b.columns)
            print(df_c.columns)
            print(f"Ready to merge.")

            # hygiene
            for k in keys:
                df_a[k] = df_a[k].astype(str).str.strip()
                df_b[k] = df_b[k].astype(str).str.strip()

            df_a = df_a.rename(
                columns={c: f"viral_{c}" for c in refseq_cols if c in df_a.columns}
            )

            df_b = df_b.rename(
                columns={c: f"bacterial_{c}" for c in refseq_cols if c in df_b.columns}
            )

            bacterial_refseq_cols = [f"bacterial_{c}" for c in refseq_cols]

            df_AA = df_a.merge(
                df_b[keys + bacterial_refseq_cols],
                on=["Sample", "SeqID", "TaxID"],
                how="left",
                validate="1:1"
            )

            print(f"AA rows: {len(df_AA)}")

            out_save_merge = f"{out_dir}/MetaFilt-Fifth-Pass.csv"
            
            df_AA["TaxID"] = df_AA["TaxID"].astype(str); df_c["TaxID"] = df_c["TaxID"].astype(str)
            df_merged = df_AA.merge(df_c[keys + bcols], on=keys, how="left")
            df_merged.to_csv(out_save_merge, index=False)
            print(f"Merged.")
            print(f"Saved to: {out_save_merge}")
            print(f"First rows (refseq columns):")
            refseq_columns_renamed = (
                [f"viral_{c}" for c in refseq_cols] +
                [f"bacterial_{c}" for c in refseq_cols]
            )
            # only return refseq_top_accession + refseq_accession_counts_json
            col_targets_substring = ["refseq_top_accession", "refseq_accession_counts_json", "accession_name"]
            refseq_columns_renamed = [c for c in refseq_columns_renamed if any(sub in c for sub in col_targets_substring)]
            subset_columns = ['name'] + keys + refseq_columns_renamed
            print(df_merged[subset_columns].head())

            # viral-only rows (has viral refseq, no bacterial)
            viral_only = df_merged[
                df_merged["viral_refseq_top_accession"].notna() &
                df_merged["bacterial_refseq_top_accession"].isna()
            ]
            save_viral_subset = f"{out_dir}/MetaFilt-Fifth-Pass-VIRAL-ONLY.csv"
            viral_only[subset_columns].to_csv(save_viral_subset, index=False)

            # bacterial-only rows (has bacterial refseq, no viral)
            bacterial_only = df_merged[
                df_merged["bacterial_refseq_top_accession"].notna() &
                df_merged["viral_refseq_top_accession"].isna()
            ]
            save_bacterial_subset = f"{out_dir}/MetaFilt-Fifth-Pass-BACTERIAL-ONLY.csv"
            bacterial_only[subset_columns].to_csv(save_bacterial_subset, index=False)

            print("\nViral-only RefSeq rows:")
            print(f"Number of viral-only rows: {len(viral_only)}")
            print(viral_only[subset_columns].head())

            print("\nBacterial-only RefSeq rows:")
            print(f"Number of bacterial-only rows: {len(bacterial_only)}")
            print(bacterial_only[subset_columns].head())
        sys.exit(0)

    for refseq in [cs_out, csb_out]:

        if refseq == cs_out:
            _type_ = 'refseq-viral'
            refseq_tsv = VIRAL_REFSEQ_TSV
        elif refseq == csb_out:
            _type_ = 'refseq-bacterial'
            refseq_tsv = BACTERIAL_REFSEQ_TSV

        # Where you want outputs
        OUT_PARENT_SUMMARY = f"{refseq}/counter-screen-{_type_}-parent-summary.csv"
        OUT_MERGED_FILTER  = f"{refseq}/MetaFilt-Third-Pass-untargeted-uviral25-2-complexity-counter-screen-{_type_}.csv"
        OUT_SAMPLE_ROLLUP  = f"{refseq}/counter-screen-{_type_}-sample-rollup.csv"

        all_parent = []
        unique_rolls = []

        for sample in samples:
            paf_glob   = f"{refseq}/{sample}/mammalian/*_taxID-*_mammalian.paf"
            stats_glob = f"{refseq}/{sample}/mammalian/*_taxID-*_mammalian_stats.csv"

            print(paf_glob)
            print(stats_glob)

            stats_df = load_parent_stats(stats_glob)

            parent_summary_sample = build_parent_summary_from_pafs(
                counter_screen_out=refseq,
                paf_glob=paf_glob,
                nm_ratio_thresh=NM_RATIO_THRESH,
                stats_df=stats_df if not stats_df.empty else None,
                tp_require_primary=TP_REQUIRE_PRIMARY,

                nm_denom=NM_DENOM,
                min_aln_len=MIN_ALN_LEN,
            )

            print(f"[INFO] {sample}: parent_summary rows = {len(parent_summary_sample)}")

            if parent_summary_sample.empty:
                print(f"[WARN] No parent summary rows for {sample} (glob/regex mismatch?)")
                continue

            # save to sample
            sample_out_path = f"{refseq}/{sample}/counter-screen-{_type_}-parent-summary.csv"
            parent_summary_sample.to_csv(sample_out_path, index=False)
            print(f"[OK] Wrote sample parent summary: {sample_out_path}")

            all_parent.append(parent_summary_sample)

            unique_roll = build_sample_unique_rollup_from_pafs(
                sample=sample,
                paf_glob=f"{refseq}/{sample}/mammalian/*_taxID-*_mammalian.paf",
                nm_ratio_thresh=NM_RATIO_THRESH,                
                tp_require_primary=TP_REQUIRE_PRIMARY,
                nm_denom=NM_DENOM,
                min_aln_len=MIN_ALN_LEN,
            )
            unique_roll.to_csv(f"{refseq}/{sample}/counter-screen-{_type_}-sample-rollup-UNIQUE.csv", index=False)
            unique_rolls.append(unique_roll)


        if not all_parent:
            print("[ERROR] No parent summaries produced across all samples.")
            sys.exit(1)

        parent_summary = pd.concat(all_parent, ignore_index=True)

        acc_map = pd.read_csv(refseq_tsv, sep="\t", header=None, names=["accession", "accession_name"], usecols=[0, 1],dtype=str,)
        # 1) Add name for the top accession
        parent_summary = parent_summary.merge(acc_map,left_on="refseq_top_accession", right_on="accession", how="left", ).drop(columns=["accession"], errors="ignore")

        # Keep blanks instead of NaN
        parent_summary["accession_name"] = parent_summary["accession_name"].fillna("")

        parent_summary.to_csv(OUT_PARENT_SUMMARY, index=False)
        print(f"[OK] Wrote parent summary: {OUT_PARENT_SUMMARY}")

        merged = merge_into_df_filter(df_filter, parent_summary)
        merged.to_csv(OUT_MERGED_FILTER, index=False)
        print(f"[OK] Wrote merged table: {OUT_MERGED_FILTER}")

        sample_roll = build_sample_rollup(parent_summary)
        sample_roll.to_csv(OUT_SAMPLE_ROLLUP, index=False)
        print(f"[OK] Wrote sample rollup: {OUT_SAMPLE_ROLLUP}")

        cat_unique_roll = pd.concat(unique_rolls, ignore_index=True)
        OUT_CAT_UNIQUE_ROLL = f"{refseq}/counter-screen-{_type_}-sample-rollup-UNIQUE.csv"
        cat_unique_roll.to_csv(OUT_CAT_UNIQUE_ROLL, index=False)
        print(f"[OK] Wrote concatenated unique sample rollup: {OUT_CAT_UNIQUE_ROLL}")




if __name__ == "__main__":
    main()

