#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Collect taxID-level alignment statistics per sample.

For each sample:
- Find BAMs matching *taxID*.0*primary*p.bam
- Extract taxID from filename
- Count total alignments
- Count primary mapped reads (tp:A:P)

Hard-coded test case included.
"""

import subprocess
import glob
import os
import json
import re
import pandas as pd
from collections import defaultdict
import subprocess


# -----------------------------
# Helpers
# -----------------------------

def run_cmd(cmd):
    return subprocess.check_output(cmd, shell=True, text=True).strip()


def count_total_alignments(bam):
    return int(run_cmd(f"samtools view {bam} | wc -l"))


def count_primary_alignments(bam):
    return int(run_cmd(f"samtools view {bam} | rg -i 'tp:A:P' | wc -l"))


def extract_taxid(bam_path):
    m = TAXID_RE.search(bam_path)
    return m.group(1) if m else None


def get_primary_readnames(bam):
    """
    Return a set of QNAMEs for primary mapped reads in a BAM.
    Matches the filtering used in calculate_bam_stats:
    - Excludes unmapped reads (-F 4)
    - Excludes duplicate-marked reads (-F 0x400)
    - Filters for primary alignments (tp:A:P)
    """
    cmd = (
        f"samtools view -F 4 -F 0x400 {bam} | "  
        f"rg -i 'tp:A:P' | "
        f"cut -f1"
    )
    out = subprocess.check_output(cmd, shell=True, text=True)
    return set(out.strip().splitlines()) if out.strip() else set()


# -----------------------------
# Hard-coded test configuration
# -----------------------------

DIRECTORY = "/mnt/usersData/ADVTIG_v3_untargeted/map-ont/"

SAMPLES = [
    "20230119_DNA_advtig-Rep1B-E50_10CFU_41_36",
    "20221222_DNA_advtig-Rep2B-E50_10CFU_41_36",    
    "20221222_DNA_advtig-Rep3B-E50_10CFU_41_36",
    "20221222_DNA_advtig-Rep4B-E50_10CFU_41_36",
    "20221222_DNA_advtig-Rep5B-E50_10CFU_41_36",
    "20221222_DNA_advtig-Rep6B-E50_10CFU_41_36",
]

PATTERN = "*taxID*.0*primary*p.bam"

TAXID_RE = re.compile(r"taxID-([0-9]+\.[0-9]+)")

# -----------------------------
# Main logic
# -----------------------------

mapping_file = f"{DIRECTORY}/sample_taxid_bams.json"

if not os.path.exists(mapping_file):

    FILTER_CSV = "/mnt/usersData/ADVTIG_v3_untargeted/MetaFilt-Fifth-Pass.csv"
    df_filt_bam_find = pd.read_csv(FILTER_CSV, usecols=["Sample", "TaxID", "SeqID"])
    df_filt_bam_find["Sample"] = df_filt_bam_find["Sample"].astype(str)
    df_filt_bam_find["TaxID"] = df_filt_bam_find["TaxID"].astype(str)
    df_filt_bam_find["SeqID"] = df_filt_bam_find["SeqID"].astype(str)

    rows = []
    sample_taxid_bams = defaultdict(lambda: defaultdict(list))

    for sample in SAMPLES:
        sample_df_filt_bam_find = df_filt_bam_find.loc[df_filt_bam_find["Sample"] == sample]
        seqids_to_find = set(sample_df_filt_bam_find["SeqID"].tolist())

        print(f"\n[INFO] Processing sample: {sample}")
        sample_dir = os.path.join(DIRECTORY, sample)
        
        bam_files = []
        for stf in seqids_to_find:
            bam_glob = f"{sample_dir}/{stf}/{PATTERN}"
            bam_files.extend(glob.glob(bam_glob))

        print(f"[INFO] {sample}: found {len(bam_files)} BAMs")

        for bam in bam_files:
            taxid = extract_taxid(bam)

            if taxid is None:
                print(f"[WARN] Could not extract taxID from {bam}")
                continue

            sample_taxid_bams[sample][taxid].append(bam)

            total = count_total_alignments(bam)
            primary = count_primary_alignments(bam)

            rows.append({
                "sample": sample,
                "taxid": taxid,
                "bam": os.path.basename(bam),
                "total_alignments": total,
                "primary_alignments": primary,
            })

    print("\n[QC] Sample → taxID → BAM counts")
    for sample, taxids in sample_taxid_bams.items():
        for taxid, bams in taxids.items():
            print(sample, taxid, len(bams))

    with open(mapping_file, "w") as fh:
        json.dump(sample_taxid_bams, fh, indent=2)
    
    print(f"[OK] Saved sample→taxID→BAM map to {mapping_file}")        

else:   
    out_csv = f"{DIRECTORY}/sample_taxid_dedup_primary_counts.csv"

    if not os.path.exists(out_csv):
        print(f"[INFO] Loading existing sample→taxID→BAM map from {mapping_file}")
        with open(mapping_file) as fh:
            sample_taxid_bams = json.load(fh)

        rows = []

        for sample, taxids in sample_taxid_bams.items():
            for taxid, bam_list in taxids.items():

                all_reads = set()

                for bam in bam_list:
                    all_reads |= get_primary_readnames(bam)

                rows.append({
                    "sample": sample,
                    "taxid": taxid,
                    "dedup_primary_reads": len(all_reads),
                    "n_bams": len(bam_list),
                })

        df = pd.DataFrame(rows).sort_values(["sample", "taxid"])
        print(f"Number of rows: {len(df)}")

        df.to_csv(out_csv, index=False)

        print(f"[OK] Wrote {out_csv}")
    else:
        print(f"Loading csv from {out_csv}")
        df = pd.read_csv(out_csv)

        FILTER_CSV = "/mnt/usersData/ADVTIG_v3_untargeted/MetaFilt-Fifth-Pass.csv"
        df_filt = pd.read_csv(FILTER_CSV, usecols=["Sample", "TaxID"])

        filt_OUT_CSV = "/mnt/usersData/ADVTIG_v3_untargeted/Deduplicated-Read-Counts-MetaFilt-Fourth-Pass-FILTERED.csv"

        df = df.rename(columns={
            "sample": "Sample",
            "taxid": "TaxID"
        })

        print(df.columns)
        print(df_filt.columns)

        # Force merge keys to string on BOTH sides
        df["Sample"] = df["Sample"].astype(str)
        df["TaxID"]  = df["TaxID"].astype(int)
        df["TaxID"]  = df["TaxID"].astype(str)

        df_filt["Sample"] = df_filt["Sample"].astype(str)
        df_filt["TaxID"]  = df_filt["TaxID"].astype(str)

        print(df)
        print(df_filt)

        df_out = df.merge(
            df_filt[["Sample", "TaxID"]],
            on=["Sample", "TaxID"],
            how="inner"
        )

        print(f"[QC] Input rows:   {len(df)}")
        print(f"[QC] Filter pairs: {len(df_filt)}")
        print(f"[QC] Output rows:  {len(df_out)}")

        assert len(df_out) > 0, "ERROR: filtering removed all rows"
        df_out = df_out.drop_duplicates(subset=["Sample", "TaxID"])
        df_out.to_csv(filt_OUT_CSV, index=False)
        print(f"[OK] Wrote filtered csv to \n{filt_OUT_CSV}")
        print(f"Number of rows: {len(df_out)}")
        