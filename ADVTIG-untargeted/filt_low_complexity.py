#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# time ./ADVTIG-untargeted/filt_low_complexity.py

import pandas as pd
import numpy as np
import argparse
import glob
import json
import re
import os
import pysam
from multiprocessing import Pool, cpu_count

# -------------------------------------------------------
#   Complexity Subfunctions
# -------------------------------------------------------
def unique_kmer_fraction_simple(seq, k=4):
    """Lightweight unique kmer fraction."""
    L = len(seq)
    if L < k:
        return 0.0
    uniq = {seq[i:i+k] for i in range(L - k + 1)}
    return len(uniq) / (L - k + 1)


def classify_complexity(seq, k):
    """
    Compute unique_kmer_fraction and assign a band.
    
    Bands:
      0 = high complexity          (>=0.30)
      1 = moderate                 (0.20–0.30)
      2 = repeat-rich              (0.10–0.20)
      3 = low complexity           (0.05–0.10)
      4 = severe low complexity    (<0.05)
    """
    kf = unique_kmer_fraction_simple(seq, k)

    if kf < 0.05:
        band = 4
    elif kf < 0.10:
        band = 3
    elif kf < 0.20:
        band = 2
    elif kf < 0.30:
        band = 1
    else:
        band = 0

    return kf, band


# -------------------------------------------------------
#   SAM Processing Subfunction
# -------------------------------------------------------
def process_sample_seqid(save_dir, sample, seqid):
    """
    Locate SAM file, extract low-complexity reads,
    save results inside mapped folder.
    """

    kmer_size = 4

    mapped_folder = f"{save_dir}/map-ont/{sample}/{seqid}/"
    bam_glob = glob.glob(f"{mapped_folder}*taxID*.0-uviral25-2-map-ont-alignment-primary-map.bam")

    if not bam_glob:
        print(f"[WARN] No BAM file for {sample}/{seqid}")
        return

    bam_path = bam_glob[0]
    print(f"[INFO] Processing BAM: {bam_path}")

    total_reads = 0
    low_reads = []

    # band counters
    band_counts = {0:0, 1:0, 2:0, 3:0, 4:0}

    bam = pysam.AlignmentFile(bam_path, "rb")

    for aln in bam.fetch(until_eof=True):
        if aln.is_secondary or aln.is_supplementary:
            continue

        seq = aln.query_sequence
        if seq is None:
            continue

        total_reads += 1

        # compute KF + band
        kf, band = classify_complexity(seq, kmer_size)
        band_counts[band] += 1

    bam.close()

    low_unique = sorted(set(low_reads))

    # save IDs
    outdir = f"{save_dir}/map-ont/{sample}/complexity_filt/"
    os.makedirs(outdir, exist_ok=True)

    outpath = f"{outdir}/{seqid}_low_complexity_reads.txt"
    with open(outpath, "w") as f:
        f.writelines(r + "\n" for r in low_unique)

    print(f"[RESULT] {total_reads} reads processed "
        f"(complexity bands: {band_counts}) "
        f"FROM SAMPLE:{sample}; SEQID:{seqid}")
    print(f"[INFO] Saved band summary to → {outpath}")

    summary = {
        "Sample": sample,
        "SeqID": seqid,
        "total_reads": total_reads,

        "kmer_size": kmer_size,
        # "low_complexity_threshold": threshold,

        "complexity_band_counts": {
            "0_high": band_counts[0],
            "1_moderate": band_counts[1],
            "2_repeat_rich": band_counts[2],
            "3_low": band_counts[3],
            "4_severe": band_counts[4]
        }
    }

    with open(f"{outdir}/{seqid}_summary.json", "w") as f:
        json.dump(summary, f)


def _process_one_bam_worker(args):
    save_dir, sample, seqid = args
    process_sample_seqid(save_dir, sample, seqid)


def parse_arguments():
    """
    Parses command line arguments.

    Returns:
    argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Identify and Extract Reads Script.")
    parser.add_argument('-o', '--save-dir', type=str, required=True, help='Path to the save directory')

    return parser.parse_args()


def main(save_dir):
    third_pass_extra_filt = f"{save_dir}/analysis/MetaFilt-Third-Pass-untargeted-uviral25-2_mapping_stats_labelled.csv"
    df_main = pd.read_csv(third_pass_extra_filt)
    
    unique_pairs = set(map(tuple, df_main[["Sample", "SeqID"]].to_numpy()))

    print(f"Number of unique pairs: {len(unique_pairs)}")

    # Prepare MP job list
    BAM_tasks = [(save_dir, sample, seqid) for sample, seqid in unique_pairs]

    # Choose number of processes (6 or auto)
    NPROC = min(10, cpu_count())

    print(f"[INFO] Launching multiprocessing with {NPROC} workers…")
    
    save_summary = f"{save_dir}/map-ont/summary.csv"
    # if not os.path.exists(save_summary):

    with Pool(processes=NPROC) as pool:
        pool.map(_process_one_bam_worker, BAM_tasks)

    summaries = []
    for path in glob.glob(f"{save_dir}/map-ont/2*/complexity_filt/*summary.json"):
        entry = json.load(open(path))

        # flatten complexity_band_counts into top-level keys
        bands = entry.pop("complexity_band_counts", {})
        for k, v in bands.items():
            entry[f"complexity_{k}"] = v

        summaries.append(entry)

    df_summary = pd.DataFrame(summaries)

    print(f"Saving to: {save_summary}")
    df_summary.to_csv(save_summary, index=False)
    # else:
    #     df_summary = pd.read_csv(save_summary)

    # Merge on both Sample and SeqID
    merged = pd.merge(
        df_main,
        df_summary,
        on=["Sample", "SeqID"],
        how="left",             # preserve all rows from the mapping stats
        suffixes=("", "_complex")
    )

    # ---------------------------------------------------------
    # 🧮 Add fractional complexity columns before saving
    # ---------------------------------------------------------
    if "total_reads" in merged.columns:
        for band in ["0_high", "1_moderate", "2_repeat_rich", "3_low", "4_severe"]:
            count_col = f"complexity_{band}"
            frac_col = f"{count_col}_frac"

            if count_col in merged.columns:
                merged[frac_col] = merged[count_col] / merged["total_reads"]
            else:
                merged[frac_col] = float("nan")

        # combined low-complexity fraction (useful diagnostic)
        merged["low_complexity_frac"] = (
            merged["complexity_3_low_frac"].fillna(0)
            + merged["complexity_4_severe_frac"].fillna(0)
        )
    else:
        print("[WARN] Column 'total_reads' not found — skipping fraction computation.")
    
    merged["HighComplexityGreater"] = merged["complexity_0_high_frac"] > merged["low_complexity_frac"]
    print(f"[INFO] Merged dataframe: {merged.shape[0]} rows, {merged.shape[1]} cols")
    merged_out = f"{save_dir}/analysis/MetaFilt-Third-Pass-untargeted-uviral25-2-complexity.csv"

    # Save merged CSV
    merged.to_csv(merged_out, index=False)
    print(f"[INFO] Saved merged file → {merged_out}")


if __name__ == "__main__":

    args = parse_arguments()
    save_dir = args.save_dir
    main(save_dir)