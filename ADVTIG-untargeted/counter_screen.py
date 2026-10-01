#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Align a single FASTA/FASTQ of reads against two hard-coded minimap2 indexes:
  1) ADV5.mmi
  2) all_mammalian.mmi

Outputs separate PAF files for each panel.
"""

import argparse
import subprocess
import os
import glob
import sys
import pandas as pd
from multiprocessing import Pool, cpu_count

# Hard-coded MMIs (adjust paths if they live elsewhere)
ADV_MMI = "/mnt/usersData/ADVTIG_fa/ADV5.mmi"
MAMMALIAN_MMI = "/mnt/usersData/ADVTIG_fa/all_mammalian.mmi"
LAMBDA_MMI = "/home/james/SequencingData/Centrifuge_libraries/ecoli-lambda-phage.NC_001416.1.mmi"

INDEX_MAMMALIAN = "/mnt/usersData/ADVTIG_fa/mammalian_ref_species2.tsv"

MAMMALIAN_ALL_MMI = "/mnt/usersData/ADVTIG_fa/all_mammalian_targets.mmi"
INDEX_MAMMALIAN_ALL = "/mnt/usersData/ADVTIG_fa/mammalian_ref_species2.tsv"

# with PERVs + LAMBDA
MAMMALIAN_ALL_MMI2 = "/mnt/usersData/ADVTIG_fa/all_mammalian_targets2.mmi"
INDEX_MAMMALIAN_ALL2 = "/mnt/usersData/ADVTIG_fa/mammalian_ref_species2.tsv"

VIRAL_REFSEQ_MMI = "/mnt/usersData/ADVTIG_fa/viral-input-sequences.mmi"
VIRAL_REFSEQ_TSV = "/mnt/usersData/ADVTIG_fa/viral-input-sequences.tsv"

BACTERIAL_REFSEQ_MMI="/mnt/usersData/ADVTIG_fa/bacterial-input-sequences.mmi"
BACTERIAL_REFSEQ_TSV="/mnt/usersData/ADVTIG_fa/bacterial-input-sequences.tsv"

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Align reads against ADV5 and all_mammalian minimap2 panels."    )
    parser.add_argument("-c", "--concatenate", action="store_true", help="Concatenate all mammalian reads per sample before alignment.")
    parser.add_argument("-o", "--out-dir", type=str, required=True, help="Directory to save alignment outputs.")
    parser.add_argument("--asm5", action="store_true", help="Use -ax asm5 instead of -ax map-ont for both panels.")    
    parser.add_argument("-s", "--samples-config", required=True, help="Sample configuration CSV (same as glue.py)")
    parser.add_argument("-p", "--processes", type=str, default=4, help="Number of parallel processes to use.")
    parser.add_argument("--panel", type=str, choices=["split", "new_combined", "old_combined", "refseq", "bacterial_refseq"], default="split", help="Whether to align to panels separately (split) or combined (combined).")
    parser.add_argument("-r", "--refseq-split",action="store_true",help="(Only for --panel refseq) align each (SeqID,TaxID) subset separately instead of sample-level concat.")
    parser.add_argument("-e", "--edit-distance-threshold", type=float, help="Edit distance threshold for filtering alignments.")

    return parser.parse_args()


def concat_sample_fastqs_for_panel(sample, seqIDs_taxIDs, out_dir, counter_screen_out, label):
    """
    Concatenate all FASTQs for a sample into one file for alignment.
    Returns: (concat_fastq_path, read_id_mapping_file)
    """
    sample_panel_dir = f"{counter_screen_out}/{sample}/{label}"
    os.makedirs(sample_panel_dir, exist_ok=True)
    
    concat_fastq = f"{sample_panel_dir}/{sample}_concat.fastq"
    read_mapping_file = f"{sample_panel_dir}/{sample}_read_mapping.csv"
    
    # Track which reads belong to which seqid/taxid
    read_mappings = []
    
    with open(concat_fastq, "w") as out_fq:
        for seqid, taxid in seqIDs_taxIDs:
            seqid = str(seqid)
            taxid = str(taxid)
            if seqid == "nan" or taxid == "nan":
                continue
            
            primary_mapped_bams = f"{out_dir}/map-ont/{sample}/{seqid}/*taxID-{taxid}.*-uviral25-2-map-ont-alignment-primary-map.bam"
            bam_glob = glob.glob(primary_mapped_bams)
            if len(bam_glob) == 0:
                continue
            if len(bam_glob) != 1:
                continue
            
            bam_file = bam_glob[0]
            
            # Convert BAM to FASTQ and track read IDs
            cmd = ["samtools", "fastq", "-F", "4", "-F", "0x400", bam_file]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            # Parse FASTQ and write to concat file while tracking read IDs
            lines = result.stdout.strip().split("\n")
            for i in range(0, len(lines), 4):
                if i + 3 < len(lines):
                    read_id = lines[i].split()[0][1:]  # Remove @ prefix
                    out_fq.write(f"{lines[i]}\n{lines[i+1]}\n{lines[i+2]}\n{lines[i+3]}\n")
                    read_mappings.append([read_id, seqid, taxid])
    
    # Save read mapping
    df_mapping = pd.DataFrame(read_mappings, columns=["read_id", "seqid", "taxid"])
    df_mapping.to_csv(read_mapping_file, index=False)
    
    print(f"Concatenated {len(read_mappings)} reads for {sample}/{label}")
    
    return concat_fastq, read_mapping_file


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


def bam_to_fastq(bam_file, fastq_file):
    print(f"Converting {bam_file} to {fastq_file}")
    # -F 4: exclude unmapped reads
    # -F 0x400: exclude secondary alignments
    cmd = ["samtools", "fastq", "-F", "4", "-F", "0x400", bam_file]
    with open(fastq_file, "w") as f:
        subprocess.run(cmd, stdout=f, check=True)
    return fastq_file


def align_to_minimap2(full_fastq: str, reference: str, asm5: bool, paf_aligned: str):
    print(f"Aligning {full_fastq} to {reference}")
    if asm5:
        minimap2 = subprocess.run(
            ["minimap2", "-cx", "asm5", reference, full_fastq],
            capture_output=True
        )
    else:
        minimap2 = subprocess.run(
            ["minimap2", "-cx", "map-ont", reference, full_fastq],
            capture_output=True
        )

    returncode = minimap2.returncode
    if returncode != 0:
        print("No alignment found (minimap2 return code:", returncode, ")")

    with open(paf_aligned, "w") as f:
        f.write(str(minimap2.stdout, "utf-8"))
        print(f"Size of PAF file: {os.stat(paf_aligned).st_size / 1e6:.0f} MB")


def run_minimap2_task(task):
    """
    Worker function for multiprocessing.
    task = [label, mmi, sample, bam_file, taxid, seqid, counter_screen_out, asm5]
    """
    label, mmi, edit_distance_threshold, sample, bam_file, taxid, seqid, counter_screen_out, asm5 = task
    
    # Create output directory for this sample/panel
    sample_panel_dir = f"{counter_screen_out}/{sample}/{label}"
    os.makedirs(sample_panel_dir, exist_ok=True)
    
    # Extract reads from BAM to FASTQ
    fastq_file = f"{sample_panel_dir}/{seqid}_taxID-{taxid}.fastq"
    bam_to_fastq(bam_file, fastq_file)
    
    # Check if FASTQ has content
    if os.stat(fastq_file).st_size == 0:
        print(f"[WARNING] Empty FASTQ for {sample}, {seqid}, {taxid}")
        return None
    
    # Count total reads
    total_reads = 0
    with open(fastq_file, "r") as f:
        for line in f:
            total_reads += 1
    total_reads = total_reads // 4
    
    # Run minimap2 to generate PAF
    paf_file = f"{sample_panel_dir}/{seqid}_taxID-{taxid}_{label}.paf"
    stats_file = f"{sample_panel_dir}/{seqid}_taxID-{taxid}_{label}_stats.csv"
    # if os.path.exists(paf_file):
    #     return

    align_to_minimap2(fastq_file, mmi, asm5, paf_file)
    
    # Parse PAF and count passing reads
    passing_reads = set()
    
    if os.path.exists(paf_file) and os.stat(paf_file).st_size > 0:
        with open(paf_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 12:
                    continue
                
                read_id = fields[0]
                query_len = int(fields[1])
                
                nm_tag = None
                for field in fields[12:]:
                    if field.startswith("NM:i:"):
                        nm_tag = int(field.split(":")[-1])
                        break
                
                if nm_tag is None:
                    continue
                
                edit_distance_ratio = nm_tag / query_len if query_len > 0 else 1.0
                
                if edit_distance_ratio <= edit_distance_threshold:
                    passing_reads.add(read_id)
    
    passing_count = len(passing_reads)
    has_hits = passing_count > 0
    passing_pct = round((passing_count / total_reads) * 100, 2) if total_reads > 0 else 0.0
    
    # Save stats
    with open(stats_file, "w") as f:
        f.write("sample,panel,seqid,taxid,has_hits,passing_reads,passing_pct\n")
        f.write(f"{sample},{label},{seqid},{taxid},{has_hits},{passing_count},{passing_pct}\n")
    
    print(f"Completed: {sample} | {label} | {seqid} | {taxid} | has_hits: {has_hits} | {passing_count}/{total_reads} ({passing_pct}%)")
    
    # Clean up intermediate files
    if os.path.exists(fastq_file):
        os.remove(fastq_file)


def fill_nan_flags_and_counts(df):
    """
    Fill NaN values for:
      - boolean-like columns (panels + species flags) -> False
      - numeric count/percentage columns (*_reads, *_pct, mammal_total_*) -> 0
    Modifies df in-place and also returns it.
    """
    # Boolean-like columns
    bool_cols = [
        col for col in df.columns
        if df[col].dtype == bool
           or set(df[col].dropna().unique()) <= {True, False}
    ]
    for col in bool_cols:
        df[col] = df[col].fillna(False)

    # Numeric count/percentage columns
    num_cols = [
        col for col in df.columns
        if (
            ("_reads" in col)
            or ("_pct" in col)
            or col in ["mammal_total_reads", "mammal_total_pct"]
        )
        and pd.api.types.is_numeric_dtype(df[col])
    ]
    for col in num_cols:
        df[col] = df[col].fillna(0)

    return df


def build_passing_reads_by_species(paf_file, ref_species_map, edit_distance_threshold):
    """
    Parse a single mammalian concat PAF and return:
      passing_reads_by_species: dict[read_id] -> set({species_1, species_2, ...})
    """
    passing_reads_by_species = {}

    if not os.path.exists(paf_file) or os.stat(paf_file).st_size == 0:
        return passing_reads_by_species

    with open(paf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 12:
                continue

            read_id = fields[0]
            query_len = int(fields[1])
            ref_name = fields[5]

            # NM tag
            nm_tag = None
            for field in fields[12:]:
                if field.startswith("NM:i:"):
                    nm_tag = int(field.split(":")[-1])
                    break
            if nm_tag is None:
                continue

            if query_len <= 0:
                continue

            edit_distance_ratio = nm_tag / float(query_len)
            if edit_distance_ratio > edit_distance_threshold:
                continue

            species = ref_species_map.get(ref_name, "Unknown")
            if species == "Unknown":
                continue

            if read_id not in passing_reads_by_species:
                passing_reads_by_species[read_id] = set()
            passing_reads_by_species[read_id].add(species)

    return passing_reads_by_species


def update_stats_files_with_species(sample_dir, read_mapping_df, passing_reads_by_species, ref_species_map, allowed_species=None):
    """
    For all *_mammalian_stats.csv in sample_dir, add per-species
    bool / reads / pct columns based on passing_reads_by_species.

    If allowed_species is not None, only those species names will get columns.
    Otherwise, all species in ref_species_map (except 'Unknown') are used.
    """
    if allowed_species is None:
        # All species labels from the TSV (mammals + viruses)
        all_species = sorted(set(ref_species_map.values()) - {"Unknown"})
    else:
        # Use only the species the caller cares about
        all_species = sorted(set(allowed_species) - {"Unknown"})

    stats_pattern = f"{sample_dir}/*_taxID-*_mammalian_stats.csv"
    stats_files = glob.glob(stats_pattern)

    for stats_file in stats_files:
        df_stats = pd.read_csv(stats_file)

        # Skip if we've already added *_reads columns once
        if any(col.endswith("_reads") for col in df_stats.columns):
            # print(f"Skipping {stats_file} - species count columns already exist")
            continue

        for idx, row in df_stats.iterrows():
            seqid = row["seqid"]
            taxid = row["taxid"]

            group = read_mapping_df[
                (read_mapping_df["seqid"] == seqid) &
                (read_mapping_df["taxid"] == taxid)
            ]
            group_read_ids = set(group["read_id"].values)
            total_reads = len(group_read_ids)

            # Initialise per-row structures
            species_hits  = {sp: False for sp in all_species}
            species_counts = {sp: 0 for sp in all_species}
            species_pcts   = {sp: 0.0 for sp in all_species}

            # Count reads per species
            for read_id in group_read_ids:
                if read_id in passing_reads_by_species:
                    for sp in passing_reads_by_species[read_id]:
                        if sp in species_hits:
                            species_hits[sp] = True
                            species_counts[sp] += 1

            # Compute percentages
            if total_reads > 0:
                for sp in all_species:
                    species_pcts[sp] = round(
                        species_counts[sp] / float(total_reads) * 100.0,
                        2
                    )

            # Write into df_stats
            for sp in all_species:
                df_stats.at[idx, sp] = species_hits[sp]
                df_stats.at[idx, f"{sp}_reads"] = species_counts[sp]
                df_stats.at[idx, f"{sp}_pct"] = species_pcts[sp]

        df_stats.to_csv(stats_file, index=False)
        print(f"Updated {stats_file}")


def add_species_columns_to_mammalian_stats(counter_screen_out, ref_species_map, allowed_species, edit_distance_threshold):
    """
    Wrapper: for each sample's mammalian concat PAF, build passing_reads_by_species
    and update its *_mammalian_stats.csv files with species columns.
    """
    paf_pattern = f"{counter_screen_out}/*/mammalian/*_concat.paf"
    paf_files = glob.glob(paf_pattern)

    for paf_file in paf_files:
        sample_dir = os.path.dirname(paf_file)
        sample = os.path.basename(os.path.dirname(sample_dir))

        read_mapping_file = paf_file.replace("_concat.paf", "_read_mapping.csv")
        if not os.path.exists(read_mapping_file):
            print(f"[WARN] Missing read_mapping for {sample}: {read_mapping_file}")
            continue

        read_mapping_df = pd.read_csv(read_mapping_file)

        passing_reads_by_species = build_passing_reads_by_species(
            paf_file,
            ref_species_map,
            edit_distance_threshold
        )

        update_stats_files_with_species(
            sample_dir,
            read_mapping_df,
            passing_reads_by_species,
            ref_species_map,
            allowed_species=allowed_species,
        )


def run_minimap2_concat_task(task):
    """
    Worker function for multiprocessing concat alignments.
    task = ["mammalian_concat", mmi, sample, concat_fastq, read_mapping_file, counter_screen_out, asm5]
    """
    label, mmi, edit_distance_threshold, sample, concat_fastq, read_mapping_file, counter_screen_out, asm5 = task
    
    sample_panel_dir = f"{counter_screen_out}/{sample}/mammalian"
    os.makedirs(sample_panel_dir, exist_ok=True)
    
    # Check if concat FASTQ has content
    if not os.path.exists(concat_fastq) or os.stat(concat_fastq).st_size == 0:
        print(f"[WARNING] Empty concat FASTQ for {sample}")
        return None
    
    # Run minimap2 on concatenated FASTQ
    paf_file = f"{sample_panel_dir}/{sample}_concat.paf"
    csv_file = paf_file.replace("_concat.paf", "_read_mapping.csv")
    if os.path.exists(paf_file) and os.path.exists(csv_file):
        return  # Skip if already done

    align_to_minimap2(concat_fastq, mmi, asm5, paf_file)
    
    # Load read mapping
    read_mapping_df = pd.read_csv(read_mapping_file)
    
    # Parse PAF and collect passing reads by read_id
    passing_reads = set()

    
    if os.path.exists(paf_file) and os.stat(paf_file).st_size > 0:
        with open(paf_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 12:
                    continue
                
                read_id = fields[0]
                query_len = int(fields[1])
                
                # Find NM tag
                nm_tag = None
                for field in fields[12:]:
                    if field.startswith("NM:i:"):
                        nm_tag = int(field.split(":")[-1])
                        break
                
                if nm_tag is None:
                    continue
                
                # Calculate edit distance ratio
                edit_distance_ratio = nm_tag / query_len if query_len > 0 else 1.0
                
                # Check if passes threshold
                if edit_distance_ratio <= edit_distance_threshold:
                    passing_reads.add(read_id)
    
    print(f"  Found {len(passing_reads)} passing reads for {sample}/mammalian")
    
    # Group by seqid/taxid and check if any reads passed
    seqid_taxid_groups = read_mapping_df.groupby(["seqid", "taxid"])
    
    for (seqid, taxid), group in seqid_taxid_groups:
        # Check if any reads from this seqid/taxid passed the filter
        group_read_ids = set(group["read_id"].values)
        has_hits = bool(group_read_ids & passing_reads)
        
        # Save stats
        stats_file = f"{sample_panel_dir}/{seqid}_taxID-{taxid}_mammalian_stats.csv"
        with open(stats_file, "w") as f:
            f.write("sample,panel,seqid,taxid,has_hits\n")
            f.write(f"{sample},mammalian,{seqid},{taxid},{has_hits}\n")
        
        print(f"Completed: {sample} | mammalian | {seqid} | {taxid} | has_hits: {has_hits}")
    
    # Clean up concat files
    if os.path.exists(concat_fastq):
        os.remove(concat_fastq)


def summarize_refseq_paf(paf_file, ref_species_map, sample_name, out_csv, edit_distance_threshold):
    """
    Summarise a concat PAF for RefSeq viral panel.

    For each reference accession (field 5 in PAF), count how many unique reads
    pass the edit-distance filter, and attach the human-readable name from
    ref_species_map (accession -> species/name).

    Output columns:
      sample, accession, name, passing_reads, passing_pct
    """
    if not os.path.exists(paf_file) or os.stat(paf_file).st_size == 0:
        print(f"[WARN] Empty or missing PAF: {paf_file}")
        return

    accession_to_reads = {}   # accession -> set(read_ids)
    all_reads = set()         # all unique query read IDs seen (for percentage)

    with open(paf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            read_id = fields[0]
            try:
                query_len = int(fields[1])
            except ValueError:
                continue

            ref_name = fields[5]   # accession in the reference index
            all_reads.add(read_id)

            if query_len <= 0:
                continue

            # Find NM:i tag
            nm_tag = None
            for field in fields[12:]:
                if field.startswith("NM:i:"):
                    try:
                        nm_tag = int(field.split(":")[-1])
                    except ValueError:
                        nm_tag = None
                    break

            if nm_tag is None:
                continue

            edit_ratio = nm_tag / float(query_len)
            if edit_ratio > edit_distance_threshold:
                continue

            # Keep one "vote" per (read, accession)
            if ref_name not in accession_to_reads:
                accession_to_reads[ref_name] = set()
            accession_to_reads[ref_name].add(read_id)

    rows = []
    total_reads = len(all_reads) if all_reads else 0

    for accession, reads in accession_to_reads.items():
        # Strip any leading ">" in keys of ref_species_map, robustly
        # (in case TSV was built that way).
        acc_key = accession.lstrip(">")
        name = ref_species_map.get(acc_key, "Unknown")
        n_reads = len(reads)
        pct = (n_reads / total_reads * 100.0) if total_reads > 0 else 0.0
        rows.append([sample_name, accession, name, n_reads, round(pct, 2)])

    if not rows:
        print(f"[INFO] No passing alignments in {paf_file}")
        return

    df_sum = pd.DataFrame(
        rows,
        columns=["sample", "accession", "name", "passing_reads", "passing_pct"],
    ).sort_values("passing_reads", ascending=False)

    # need to add a step here to add the taxID based on accession

    df_sum.to_csv(out_csv, index=False)
    print(f"[REFSEQ] Wrote summary for {sample_name} to {out_csv}")
    return df_sum


def main():
    args = parse_arguments()
    samples_config = args.samples_config
    out_dir = args.out_dir
    asm5 = args.asm5
    processes = args.processes
    concatenate = args.concatenate
    panel = args.panel
    refseq_split = args.refseq_split

    edit_distance_threshold = args.edit_distance_threshold

    print(f"\nRunning counter screen with panel: {panel}, concatenate: {concatenate}, refseq_split: {refseq_split}, edit_distance_threshold: {edit_distance_threshold}\n")

    merged_out = f"{out_dir}/analysis/MetaFilt-Third-Pass-untargeted-uviral25-2-complexity.csv"
    allowed_species = None

    if panel == "new_combined":
        counter_screen_out = f"{out_dir}/counter_screen_comb3/"
        save_summary = merged_out.replace(".csv", f"-counter-screened-summary-combined3.csv")
        panels = [
            ("mammalian", MAMMALIAN_ALL_MMI2),
        ]
        print(INDEX_MAMMALIAN_ALL2)
        ref_species_df = pd.read_csv(INDEX_MAMMALIAN_ALL2, sep="\t", names=["accession", "species"])
    elif panel == "old_combined":
        counter_screen_out = f"{out_dir}/counter_screen_comb/"
        save_summary = merged_out.replace(".csv", f"-counter-screened-summary-combined.csv")
        panels = [
            ("mammalian", MAMMALIAN_ALL_MMI),
        ]
        ref_species_df = pd.read_csv(INDEX_MAMMALIAN_ALL, sep="\t", names=["accession", "species"])
    elif panel == "refseq":
        counter_screen_out = f"{out_dir}/counter_screen_refseq/"
        save_summary = merged_out.replace(".csv", f"-counter-screened-summary-refseq.csv")
        panels = [
            ("mammalian", VIRAL_REFSEQ_MMI),
        ]
        ref_species_df = pd.read_csv(VIRAL_REFSEQ_TSV, sep="\t", names=["accession", "species"])
        # Only keep mammalian + targeted viruses in the summary
        canonical_species = {
            "Cat", "Human", "Monkey", "Pig", "Tamarin",
            "RSV_A", "REOvirus", "PCV1", "NC_001510_1", "AF038600.1", "DNA_CS",
            "EBV_B95_8", "FeLV", "ADV5", "NC_001514_1",
            "FeLV_Kawakami_Theilen", "HCoV_NL63",
        }
        # restrict to species that actually exist in the TSV
        allowed_species = sorted(set(ref_species_df["species"].unique()) & canonical_species)
    elif panel == "bacterial_refseq":
        counter_screen_out = f"{out_dir}/counter_screen_bacterial_refseq/"
        save_summary = merged_out.replace(".csv", f"-counter-screened-summary-bacterial-refseq.csv")
        panels = [
            ("mammalian", BACTERIAL_REFSEQ_MMI),
        ]
        ref_species_df = pd.read_csv(BACTERIAL_REFSEQ_TSV, sep="\t", names=["accession", "species"])        
    else:
        counter_screen_out = f"{out_dir}/counter_screen/"
        save_summary = merged_out.replace(".csv", f"-counter_screened_summary.csv")
        panels = [
            ("DNA_CS",    LAMBDA_MMI),
            ("ADV5",      ADV_MMI),
            ("mammalian", MAMMALIAN_MMI),
        ]
        ref_species_df = pd.read_csv(INDEX_MAMMALIAN, sep="\t", names=["accession", "species"])

    samples, _, _ = load_samples(samples_config)
    print(f"Samples: {samples}")
    
    os.makedirs(counter_screen_out, exist_ok=True)

    df = pd.read_csv(merged_out)
    df_filter = df.loc[df["low_complexity_frac"] <= 0.2] # change from <0.2 to <=0.2
    df_filter = df_filter.loc[df_filter["Depth 1X"] >= 0.6]           # changed from 80% -> 70% ## change to 60%?
    df_filter = df_filter.loc[df_filter["Identity Percentage"] >= 60] # changed from 80% -> 70% ## change to 60%?
    df_filter = df_filter.loc[df_filter["Mapped reads"] > 1] 
    df_filter = df_filter.loc[df_filter["ED_passing_ratio"] >= 0.6] # changed from 70% -> 60%
    print(f"Total targets before filtering: {len(df)} and after filtering: {len(df_filter)}")

    all_runs = []
    cat_runs = []
    for label, mmi in panels:
        print(f"\n\nProcessing panel: {label} with MMI: {mmi}")
        for sample in samples:
            print(f"Current sample: {sample}")
            sub_df = df_filter.loc[df_filter["Sample"] == sample]
            if sub_df.empty:
                print(f"[WARNING] No data for sample {sample} in {merged_out}, skipping minimap2 alignment.")
                continue

            seqIDs_taxIDs = sub_df[["SeqID", "TaxID"]].drop_duplicates().to_numpy()
            print(f"Number of unique (SeqID, TaxID) pairs: {len(seqIDs_taxIDs)}")

            # add handling for bacterial_refseq panel
            if label == "mammalian" and concatenate and not (
                (panel == "refseq" and refseq_split) or
                (panel == "bacterial_refseq" and refseq_split)
            ):
            # if label == "mammalian" and concatenate and not (panel == "refseq" and refseq_split):
                concat_fastq, read_mapping = concat_sample_fastqs_for_panel(sample, seqIDs_taxIDs, out_dir, counter_screen_out, label)
                cat_runs.append(["mammalian_concat", mmi, edit_distance_threshold, sample, concat_fastq, read_mapping, counter_screen_out, asm5])
            else:  
                for seqid, taxid in seqIDs_taxIDs:
                    seqid = str(seqid)
                    # taxid = str(taxid)
                    taxid = str(int(float(taxid)))
                    if seqid == "nan" or taxid == "nan":
                        continue

                    primary_mapped_bams = f"{out_dir}/map-ont/{sample}/{seqid}/*taxID-{taxid}.*-uviral25-2-map-ont-alignment-primary-map.bam"
                    bam_glob = glob.glob(primary_mapped_bams)
                    if len(bam_glob) == 0:
                        print(f"[WARNING] No primary mapped BAM files found for sample {sample}, SeqID {seqid}, TaxID {taxid}.")
                        continue
                    if len(bam_glob) == 1:
                        all_runs.append([label, mmi, edit_distance_threshold, sample, bam_glob[0], taxid, seqid, counter_screen_out, asm5])
    
    print(f"\n\nAll runs collected: {len(all_runs)} --- concatenated runs: {len(cat_runs)}")

    with Pool(processes=int(processes)) as pool:
        pool.map(run_minimap2_task, all_runs)

    with Pool(processes=6) as pool:
        pool.map(run_minimap2_concat_task, cat_runs)

    print(f"Loading mammalian reference species mapping from: {INDEX_MAMMALIAN}")
    ref_species_map = dict(zip(ref_species_df["accession"].str.replace(">", ""), ref_species_df["species"]))
    print(f"Loaded {len(ref_species_map)} reference species mappings.")

    if concatenate and not panel in ["bacterial_refseq", "refseq"]:
        add_species_columns_to_mammalian_stats(counter_screen_out, ref_species_map, allowed_species, edit_distance_threshold)

    if panel == "refseq" or panel == "bacterial_refseq":
        df_filter2 = df_filter.copy()
        df_seqIDs_taxIDs = df_filter2[["SeqID", "TaxID"]].drop_duplicates()

        print(f"\n\n{panel}:\n")
        if refseq_split:
            if panel == "refseq":
                save_summary = merged_out.replace(".csv", "-counter-screened-summary-refseq-split.csv")
            elif panel == "bacterial_refseq":
                save_summary = merged_out.replace(".csv", "-counter-screened-summary-bacterial-refseq-split.csv")
            glob_analyses = glob.glob(f"{counter_screen_out}/2*/**/*taxID*stats.csv", recursive=True)

            for c in ["Sample", "SeqID", "TaxID"]:
                df_filter2[c] = df_filter2[c].astype(str)

            all_stats = []
            for stats_file in glob_analyses:
                if not os.path.exists(stats_file):
                    continue
                df_stats = pd.read_csv(stats_file)

                # Normalize types for stable merge keys
                for c in ["sample", "seqid", "taxid"]:
                    if c in df_stats.columns:
                        df_stats[c] = df_stats[c].astype(str)

                all_stats.append(df_stats)

            stats_df = pd.concat(all_stats, ignore_index=True)
            stats_df = stats_df[["sample", "seqid", "taxid", "has_hits", "passing_reads", "passing_pct"]].copy()
            if panel == "refseq":
                stats_df = stats_df.rename(
                    columns={
                        "has_hits": "refseq_has_hits",
                        "passing_reads": "refseq_passing_reads",
                        "passing_pct": "refseq_passing_pct",
                    }
                )
                cat_df_sums = df_filter2.merge(
                    stats_df,
                    left_on=["Sample", "SeqID", "TaxID"],
                    right_on=["sample", "seqid", "taxid"],
                    how="left",
                ).drop(columns=["sample", "seqid", "taxid"], errors="ignore")
                cat_df_sums["refseq_has_hits"] = cat_df_sums["refseq_has_hits"].fillna(False)
                cat_df_sums["refseq_passing_reads"] = cat_df_sums["refseq_passing_reads"].fillna(0).astype(int)
                cat_df_sums["refseq_passing_pct"] = cat_df_sums["refseq_passing_pct"].fillna(0.0)
                print("\nRefSeq individual summarisation complete.")
            elif panel == "bacterial_refseq":
                stats_df = stats_df.rename(
                    columns={
                        "has_hits": "bacterial_refseq_has_hits",
                        "passing_reads": "bacterial_refseq_passing_reads",
                        "passing_pct": "bacterial_refseq_passing_pct",
                    }
                )
                cat_df_sums = df_filter2.merge(
                    stats_df,
                    left_on=["Sample", "SeqID", "TaxID"],
                    right_on=["sample", "seqid", "taxid"],
                    how="left",
                ).drop(columns=["sample", "seqid", "taxid"], errors="ignore")
                cat_df_sums["bacterial_refseq_has_hits"] = cat_df_sums["bacterial_refseq_has_hits"].fillna(False)
                cat_df_sums["bacterial_refseq_passing_reads"] = cat_df_sums["bacterial_refseq_passing_reads"].fillna(0).astype(int)
                cat_df_sums["bacterial_refseq_passing_pct"] = cat_df_sums["bacterial_refseq_passing_pct"].fillna(0.0)
                print("\nBacterial RefSeq individual summarisation complete.")
        else:
            df_sums = []
            for sample in samples:
                paf_file = f"{counter_screen_out}/{sample}/mammalian/{sample}_concat.paf"
                if not os.path.exists(paf_file):
                    print(f"[WARN] No concat PAF for {sample}: {paf_file}")
                    continue

                if panel == "refseq":
                    print(f"Summarizing RefSeq viral panel for sample: {sample}")
                    out_csv = paf_file.replace("_concat.paf", "_refseq_summary.csv")
                elif panel == "bacterial_refseq":
                    print(f"Summarizing Bacterial RefSeq panel for sample: {sample}")
                    out_csv = paf_file.replace("_concat.paf", "_bacterial_refseq_summary.csv")

                df_sum = summarize_refseq_paf(paf_file, ref_species_map, sample, out_csv, edit_distance_threshold)
                df_sums.append(df_sum)

            print(f"\{panel} panel summarisation complete.")
            cat_df_sums = pd.concat(df_sums)

        # if no taxID column in cat_df_sums use df_seqIDs_taxIDs to add it
        print(cat_df_sums)
        if "TaxID" not in cat_df_sums.columns:
            if "taxid" in cat_df_sums.columns:
                cat_df_sums = cat_df_sums.rename(columns={"taxid": "TaxID"})
            else:
                rename_keys = {"accession":"SeqID", "sample":"Sample", "seqID":"SeqID"}
                cat_df_sums = cat_df_sums.rename(columns=rename_keys)
                print(cat_df_sums.columns)

                # safest mapping: per-sample SeqID -> TaxID
                map_df = df_filter[["Sample", "SeqID", "TaxID"]].drop_duplicates().copy()
                print(map_df.columns)

                # normalize dtypes to avoid merge type errors
                for c in ["Sample", "SeqID", "TaxID"]:
                    map_df[c] = map_df[c].astype(str)

                for c in ["Sample", "SeqID"]:
                    if c in cat_df_sums.columns:
                        cat_df_sums[c] = cat_df_sums[c].astype(str)

                cat_df_sums = cat_df_sums.merge(
                    map_df,
                    on=["Sample", "SeqID"],
                    how="left",
                )
                # still missing majority of taxIDs.
        cat_df_sums.to_csv(save_summary, index=False)
        print(f"Refseq save: {save_summary}")

        if panel == "refseq" and not refseq_split:
            sys.exit(0)
        elif panel == "bacterial_refseq" and not refseq_split:
            sys.exit(0)

    if refseq_split:
        sys.exit(0)

    glob_analyses = glob.glob(f"{counter_screen_out}/2*/**/*taxID*stats.csv")
    all_stats = []
    for stats_file in glob_analyses:
        if os.path.exists(stats_file):
            df_stats = pd.read_csv(stats_file)
            all_stats.append(df_stats)

    if all_stats:
        # Concatenate all stats
        stats_df = pd.concat(all_stats, ignore_index=True)
        
        # For mammalian, need to handle species columns separately
        mammalian_stats = stats_df[stats_df["panel"] == "mammalian"].copy()
        other_stats = stats_df[stats_df["panel"] != "mammalian"].copy()
        
        # Pivot non-mammalian (just has_hits)
        if not other_stats.empty:
            # Keep all columns, not just has_hits
            pivot_other = other_stats.pivot_table(
                index=["sample", "seqid", "taxid"],
                columns="panel",
                values=["has_hits", "passing_reads", "passing_pct"],
                aggfunc="first"
            )
            # Flatten column names: (has_hits, ADV5) -> ADV5, (passing_reads, ADV5) -> ADV5_passing_reads
            pivot_other.columns = [f"{col[1]}_{col[0]}" if col[0] != "has_hits" else col[1] for col in pivot_other.columns]
            pivot_other = pivot_other.reset_index()
        else:
            pivot_other = None
        
        # For mammalian, keep species columns
        if not mammalian_stats.empty:
            mammalian_stats = mammalian_stats.drop(columns=["panel"], errors="ignore")
            mammalian_stats = mammalian_stats.rename(columns={"has_hits": "mammalian"})
        
        # Merge back with original filtered dataframe
        final_df = df_filter.copy()
        
        if pivot_other is not None:
            final_df = final_df.merge(
                pivot_other,
                left_on=["Sample", "SeqID", "TaxID"],
                right_on=["sample", "seqid", "taxid"],
                how="left"
            )
            final_df = final_df.drop(columns=["sample", "seqid", "taxid"], errors="ignore")
        
        if not mammalian_stats.empty:
            final_df = final_df.merge(
                mammalian_stats,
                left_on=["Sample", "SeqID", "TaxID"],
                right_on=["sample", "seqid", "taxid"],
                how="left"
            )
            final_df = final_df.drop(columns=["sample", "seqid", "taxid"], errors="ignore")
        
        final_df = fill_nan_flags_and_counts(final_df)
        
        for col in ["passing_reads", "passing_pct"]:
            if col in final_df.columns:
                final_df = final_df.drop(columns=[col])

        # Only compute mammal totals if those species columns actually exist
        mammal_species = ["Cat", "Human", "Monkey", "Pig", "Tamarin"]
        mammal_read_cols = [
            f"{sp}_reads"
            for sp in mammal_species
            if f"{sp}_reads" in final_df.columns
        ]

        if mammal_read_cols:
            final_df["mammal_total_reads"] = final_df[mammal_read_cols].sum(axis=1)
            final_df["mammal_total_pct"] = (
                final_df["mammal_total_reads"] / final_df["Mapped reads"] * 100
            )

        final_df.to_csv(save_summary, index=False)
        print(f"\n\nSaved summary to: {save_summary}")
        print(f"Total rows: {len(final_df)}")
    else:
        print("\n\nNo stats files generated")

    print("\nCounter screening complete!")


if __name__ == "__main__":
    main()
