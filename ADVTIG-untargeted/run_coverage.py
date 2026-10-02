#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pysam
import pandas as pd
import argparse
import glob
import os
import multiprocessing as mp
import json
from Bio import pairwise2, SeqIO
from Bio.Align import PairwiseAligner
import sys
import shutil
import subprocess
import re
import numpy as np
from itertools import zip_longest
from multiprocessing import Pool
from typing import Iterable, Optional, Tuple, Dict, List


def check_file_integrity(bam_file):
    result = subprocess.run(["samtools", "quickcheck", bam_file], capture_output=True)
    if result.returncode != 0:
        print(f"File {bam_file} is corrupted or invalid!!!!!!!!!!!!!!!!!!!!")
    else:
        print(f"File {bam_file} is valid.")


def extract_read_from_fastq(fastq_file_path: str, read_id: str, output_path: str):
    with open(fastq_file_path, 'r') as fastq_file, open(output_path, 'w') as output_file:
        lines = iter(fastq_file)
        for line in lines:
            if line.startswith(f'@{read_id}'):
                output_file.write(line)  # Write the read ID line
                output_file.write(next(lines))  # Write the sequence line
                output_file.write(next(lines))  # Write the '+' line
                output_file.write(next(lines))  # Write the quality line
                break


def get_longest_read(bam_file_path: str) -> str:
    command = f"samtools view {bam_file_path} | awk '{{print length($10), $0}}' | sort -nr | head -1"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.stdout:
        return result.stdout.strip().split('\n')[0]
    return None


def get_longest_read_id(bam_file_path: str) -> str:
    command = f"samtools view {bam_file_path} | awk '{{print length($10), $1}}' | sort -nr | head -1"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.stdout:
        return result.stdout.strip().split('\n')[0].split()[1]  # Returning only the read ID
    return None


def calculate_coverage_and_mean(bam_file_path: str, ref_length: int):
    """
    Calculate the total covered bases and mean coverage for a BAM file.

    Args:
        bam_file_path (str): Path to the BAM file.
        ref_length (int): Total length of the reference sequence.

    Returns:
        tuple: Total covered bases and mean coverage.
    """
    if not os.path.isfile(bam_file_path):
        print(f"Error: BAM file does not exist: {bam_file_path}")
        return 0, 0

    try:
        # Calculate total covered bases using samtools depth and awk
        command = f"samtools depth -aa {bam_file_path} | awk '{{sum+=$3}} END {{print sum}}'"
        result = subprocess.run(command, shell=True, text=True, capture_output=True)

        if result.returncode != 0:
            print(f"Error running samtools depth: {result.stderr.strip()}")
            return 0, 0

        total_covered_bases = int(result.stdout.strip()) if result.stdout.strip().isdigit() else 0

        # Calculate mean coverage
        mean_coverage = total_covered_bases / ref_length if ref_length > 0 else 0

        return total_covered_bases, mean_coverage

    except Exception as e:
        print(f"An error occurred during coverage calculation: {e}")
        return 0, 0


def get_reference_length(bam_file_path: str) -> int:
    """
    Calculate the total reference length from a BAM file.
    If the header is missing, calculate it from the alignment data.
    """
    if not os.path.isfile(bam_file_path):
        print(f"Error: BAM file does not exist: {bam_file_path}")
        return 0

    # Check if the BAM file has a header
    header_check_command = f"samtools view -H {bam_file_path}"
    header_result = subprocess.run(header_check_command, shell=True, text=True, capture_output=True)
    if header_result.stdout.strip():  # If header exists
        # Extract reference lengths from the header
        command = f"samtools view -H {bam_file_path} | grep '^@SQ' | awk '{{for(i=1;i<=NF;i++){{if($i ~ /^LN:/){{len+=substr($i,4)}}}}}} END {{print len}}'"
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        if result.stdout.strip():
            return int(result.stdout.strip())
        else:
            print("Error: Unable to parse reference length from header.")
            return 0
    else:  # Header is missing, calculate from alignment data
        print(f"Header missing for {bam_file_path}. Calculating reference length from alignment data.")
        command = f"samtools view {bam_file_path} | awk '{{if (!($1 ~ /^@/)) {{print $3}}}}' | sort | uniq -c | awk '{{print $2}}'"
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        references = result.stdout.strip().split("\n")
        total_length = 0
        for ref in references:
            # Extract max position for each reference
            pos_command = f"samtools view {bam_file_path} | awk '$3 == \"{ref}\" {{print $4}}' | sort -n | tail -1"
            pos_result = subprocess.run(pos_command, shell=True, text=True, capture_output=True)
            if pos_result.stdout.strip():
                total_length += int(pos_result.stdout.strip())
        return total_length


def calculate_bam_stats(bam_file_path: str, spp: str, sample: str) -> Dict[str, float]:
    # Helper function to run shell commands
    def run_command(command: str) -> str:
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        return result.stdout.strip()

    # Get total read count (mapped and unmapped)
    total_reads_command = f"samtools view -c {bam_file_path}"
    total_reads = int(run_command(total_reads_command))

    # Get mapped read count
    mapped_reads_command = f"samtools view -c -F 4 -F 0x400 {bam_file_path}"
    mapped_reads = int(run_command(mapped_reads_command))

    # New command for total length of mapped reads
    total_mapped_bases_command = (
        f"samtools view -F 4 -F 0x400 {bam_file_path} | "
        "awk '{SUM+=length($10)} END {print SUM}'"
    )
    total_mapped_bases_output = run_command(total_mapped_bases_command)
    total_mapped_bases = int(total_mapped_bases_output) if total_mapped_bases_output else 0

    # New calculation for mean read length of mapped reads
    mean_mapped_read_length = total_mapped_bases / mapped_reads if mapped_reads > 0 else 0

    # Command to find maximum read length of mapped reads
    max_read_length_command = (
        f"samtools view -F 4 -F 0x400 {bam_file_path} | "
        "awk '{print length($10)}' | sort -n | tail -1"
    )
    max_read_length_output = run_command(max_read_length_command)
    max_read_length = int(max_read_length_output) if max_read_length_output else 0

    # Calculate % Mapped Reads
    percent_mapped = (mapped_reads / total_reads) * 100 if total_reads > 0 else 0

    # Get the sum of all base lengths (mapped and unmapped)
    total_bases_command = f"samtools view {bam_file_path} | awk '{{SUM+=length($10)}} END {{print SUM}}'"
    total_bases_output = run_command(total_bases_command)
    total_bases = int(total_bases_output) if total_bases_output.isdigit() else 0

    # Get quality filtered read count (e.g., MAPQ >= 30)
    quality_reads_command = f"samtools view -c -q 30 {bam_file_path}"
    quality_reads = int(run_command(quality_reads_command))

    # Calculate % Mapped Quality Reads
    percent_quality = (quality_reads / total_reads) * 100 if total_reads > 0 else 0

    # Output results
    stats = {
        "Species": spp,
        "Sample": sample,
        "Total reads": total_reads,
        "Mapped reads": mapped_reads,
        "Percent mapped": percent_mapped,
        "Total bases": total_bases,
        'Max read length (mapped reads)': max_read_length,
        "Total mapped bases": total_mapped_bases,
        "Mean mapped read length": mean_mapped_read_length,
        "Quality reads (MAPQ>=30)": quality_reads,
        "Percent quality mapped": percent_quality
    }

    return stats


# Convert NumPy data types to native Python types
def convert_to_serializable(obj):
    if isinstance(obj, (np.integer, np.int64)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64)):
        return float(obj)
    elif isinstance(obj, (np.ndarray, list)):
        return [convert_to_serializable(o) for o in obj]
    else:
        return obj


def convert_to_base_pairs(value):
    if value is None:
        print("\n\n\n", value)
        return None
    # Remove 'bp' from the end of the value
    value = value.rstrip('bp')

    # Check the last character of the value to find the unit
    unit = value[-1]

    if unit == 'K':
        # The value is in kilobases. Convert it to base pairs.
        return float(value[:-1]) * 10**3
    elif unit == 'M':
        # The value is in megabases. Convert it to base pairs.
        return float(value[:-1]) * 10**6
    else:
        # The value is already in base pairs.
        return float(value)


def get_covered_bases(primary_bam):
    """
    Calculate the number of covered bases using samtools depth.
    """
    if not os.path.isfile(primary_bam):
        print(f"Error: BAM file does not exist: {primary_bam}")
        return 0, 0

    # Run samtools depth and capture the output
    command = ["samtools", "depth", "-aa", primary_bam]
    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0 or result.stderr:
        print(f"Error running samtools depth: {result.stderr}")
        return 0, 0

    # Parse samtools depth output to calculate covered bases and total bases
    covered_bases = 0
    total_cov_bases = 0
    for line in result.stdout.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) == 3:
            coverage = int(parts[2])  # Depth at this position
            total_cov_bases += coverage
            if coverage > 0:
                covered_bases += 1

    print(f"Covered bases: {covered_bases}")
    print(f"Total coverage: {total_cov_bases}")

    return covered_bases, total_cov_bases


def get_total_reference_length(bam_file):
    # Use samtools view -H to get header lines, then grep to filter lines containing reference sequence lengths
    command = f"samtools view -H {bam_file} | awk '/^@SQ/{{for(i=1;i<=NF;i++) if ($i ~ /^LN:/) print substr($i, 4)}}'"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.stderr:
        print(f"Error: {result.stderr}")
        return None

    try:
        # Extract lengths directly from the command output
        lengths = [int(line) for line in result.stdout.strip().split('\n') if line]
        total_length = sum(lengths)
        # print(f"Extracted lengths: {lengths}")
        print(f"Total reference length: {total_length}")
    except ValueError as e:
        print(f"An error occurred during length extraction: {e}")
        return None

    return total_length


def cov_analysis(db_folder, ce_name: str, database: str, sample: str, save_dir: str, fq_seqID: str, seqid_df: pd.DataFrame, combined_dict: Dict, viral_taxIDs: List, vertebrate_taxIDs: List) -> None:
    print(f"\nProcessing {sample}")
    cent_out = f"{save_dir}/analysis/sample_data/{sample}/centrifuge/"
    troubleshooting = f"{cent_out}{sample}_{ce_name}_centrifuge_troubleshooting_report.tsv"
    if not os.path.exists(troubleshooting):
        print(f"Error: {troubleshooting} does not exist")
        return
    ts_df = pd.read_csv(troubleshooting, sep="\t")

    set_seqIDs = set(list(ts_df["seqID"]))

    other_viral_paths = []
    other_viral_taxIDs = []
    for k, v in combined_dict.items():
        
        # check if k in set_seqIDs
        seqID_found = False
        if k in set_seqIDs:
            seqID_found = True

        if seqID_found:
            taxIDs = ts_df.loc[ts_df["seqID"]==k]["taxID"].unique()
            if len(taxIDs) > 0:
                taxID_found = int(taxIDs[0])
                other_viral_taxIDs.append(taxID_found)
                other_viral_path = f"/mnt/usersData//NCBI_RVDB/virus_RVDBs/{db_folder}/{v}"
                other_viral_paths.append(other_viral_path)

    # update viral_taxIDs
    # first10_taxIDs = other_viral_taxIDs[:10]
    # print(first10_taxIDs)
    viral_taxIDs = list(set(viral_taxIDs + other_viral_taxIDs))

    ts_df_known = ts_df.loc[ts_df["seqID"] != "unclassified"]
    print(f"\nUnique seqIDs: {len(ts_df_known['seqID'].unique())} --- additional viral taxIDs {len(other_viral_taxIDs)}")

    # rename code col to seqID
    seqid_df.rename(columns={"code": "seqID"}, inplace=True)
    # use seqID_df to merge with ts_df_known on seqID
    # ts_df["taxID"] = ts_df["taxID"].astype(int)
    # ts_df_known["taxID"] = ts_df_known["taxID"].astype(int)

    ts_df_known = pd.merge(ts_df_known, seqid_df, on="seqID")

    # save to json
    os.makedirs(fq_seqID, exist_ok=True)

    viral_save = f"{fq_seqID}{database}.tsv"
    vertebrate_save = f"{fq_seqID}uvertebrate.tsv"

    # subset viral
    # DOES THIS NEED MODIFYING?
    # THIS IMPLIES THAT ONLY TARGETED VIRUSES ARE BEING SUBSET AND OTHERS IGNORED
    # ts_df["taxID"] = ts_df["taxID"].astype(int)
    ts_df_known["taxID"] = ts_df_known["taxID"].astype(int)
    # ts_df_known["taxID"] = pd.to_numeric(ts_df_known["taxID"], errors="coerce").astype("Int64")

    # print(ts_df_known)
    # print(ts_df_known.loc[ts_df_known["taxID"]==47929])
    # print(ts_df_known.loc[ts_df_known["taxID"]==10243])
    # print(viral_taxIDs)
    # print("47929", 47929 in viral_taxIDs)
    # print("10243", 10243 in viral_taxIDs)
    # print("28285", 28285 in viral_taxIDs)
    # print(ts_df_known["taxID"].dtype)  # probably 'object'
    # print(ts_df_known["taxID"].apply(type).value_counts())

    # if not os.path.exists(viral_save):
    viral_df = ts_df_known.loc[ts_df_known["taxID"].isin(viral_taxIDs)]
    viral_df.to_csv(viral_save, index=False)
    print(f"Saved {viral_save}")
    # else:
    #     viral_df = pd.read_csv(viral_save)

    # if not os.path.exists(vertebrate_save):
    vertebrate_df = ts_df_known.loc[ts_df_known["taxID"].isin(vertebrate_taxIDs)]
    vertebrate_df.to_csv(vertebrate_save, index=False)
    print(f"Saved {viral_save}")
    # else:
    #     vertebrate_df = pd.read_csv(vertebrate_save)

    # get unique seqIDs
    viral_seqIDs = viral_df["seqID"].unique()
    vertebrate_seqIDs = vertebrate_df["seqID"].unique()

    # print(viral_df.loc[viral_df["taxID"]==47929])
    # first10_taxIDs_df = viral_df.loc[viral_df["taxID"].isin(first10_taxIDs)]
    print(f"\n{viral_df.head()}")#\n\n{first10_taxIDs_df}\n")
    # sys.exit(0)

    print(f"Viral sequence IDs: {len(viral_seqIDs)}")
    print(f"Vertebrate sequence IDs: {len(vertebrate_seqIDs)}")
    return viral_df, vertebrate_df, viral_seqIDs, vertebrate_seqIDs


def grouper(iterable: Iterable, n: int, fillvalue: Optional[str] = None) -> Iterable[Tuple]:
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def align_to_minimap2(full_fastq: str, reference: str, asm5: bool, sam_aligned):
    print(f"Aligning {full_fastq} to {reference}")
    if asm5:
        minimap2 = subprocess.run(["minimap2", "-ax", "asm5", reference, full_fastq], capture_output=True)
        print("minimap2", "-ax", "asm5", reference, full_fastq)
    else:
        minimap2 = subprocess.run(["minimap2", "-ax", "map-ont", reference, full_fastq], capture_output=True)
    returncode = minimap2.returncode
    if returncode != 0:
        print("No alignment found")

    with open(sam_aligned, "w") as f:
        f.write(str(minimap2.stdout, "utf-8"))
        print(f"Size of sam file: {os.stat(sam_aligned).st_size/1e6:.0f} MB")


def sam_to_bam(sam_aligned: str, bam_aligned: str):
    print(f"Converting {sam_aligned} to {bam_aligned}")
    S2B = subprocess.run(["samtools", "view", "-b", sam_aligned], capture_output=True)
    with open(bam_aligned, "wb") as ff:
        ff.write(S2B.stdout)


def sort_bam(bam_aligned: str, reference: str, sorted_path: str, mpileup_file: str):
    print(f"Sorting {bam_aligned} to {sorted_path}")
    subprocess.run(["samtools", "sort", "-o", sorted_path, bam_aligned], capture_output=True)
    print(f"Running mpileup on {sorted_path}")
    mpile_out = subprocess.run(["samtools", "mpileup", "-f", reference, sorted_path], capture_output=True)
    mPIPE = str(mpile_out.stdout, "utf-8")
    print(f"Saving mpileup to {mpileup_file}")
    with open(mpileup_file, "w") as g:
        g.write(mPIPE)
    print(f"{os.stat(mpileup_file).st_size/1e6:.0f} MB for {mpileup_file}")


def flagstat_run_and_save(sorted_path: str, flag_file: str):
    try:
        print(f"Running samtools flagstat on: {sorted_path}")
        with open(flag_file, 'w') as outfile:
            result = subprocess.run(["samtools", "flagstat", sorted_path], text=True, stdout=outfile, stderr=subprocess.PIPE, check=True)
        print(f"Result saved to {flag_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running samtools flagstat on {sorted_path}: {e.stderr}")


def primary_map(sorted_path, prim_map) -> None:
    print(f"Extracting primary mapped reads from {sorted_path}")
    mapped = subprocess.run(["samtools", "view", "-b", "-F", "0x400", "-F", "0x800", "-F", "0x100", sorted_path], capture_output=True)
    with open(prim_map, "wb") as fff:
        fff.write(mapped.stdout)
    print(f"{os.stat(prim_map).st_size / 1e3:.0f} KB for {prim_map}")


def coverage(prim_map: str, txt_map: str) -> None:
    coverage_map = subprocess.run(["samtools", "coverage", "-H", "-w", "40", prim_map], capture_output=True)
    PIPE = str(coverage_map.stdout, "utf-8")
    with open(txt_map, "w") as ffff:
        ffff.write(PIPE)


def generate_fq_seqID_files(fq_name, fq_seqID_unique, sample, directory, tsk, spp_name, readIDs):
    fq_filename = glob.glob(fq_name)
    if len(fq_filename) == 0:
        print(f"No fastq file found for {sample}")
        return
    else:
        fq_filename = fq_filename[0]
        # print(f"Found {fq_filename}")

    pattern = re.compile(r'@([\w-]+)')
    fq_seqID_dir = f"{directory}/analysis/sample_data/{sample}/{fq_seqID_unique}"
    os.makedirs(fq_seqID_dir, exist_ok=True)

    spp_name = spp_name.replace(" ", "_")
    fq_seqID_file = f"{fq_seqID_dir}/{tsk}_{spp_name}.fastq"
    print(f"Generating {fq_seqID_file}")
    if not os.path.isfile(fq_seqID_file):
        print(f"Writing to {fq_seqID_file} --- number of readIDs: {len(readIDs)}")
        with open(fq_seqID_file, "w") as outfile:
            with open(fq_filename, "r") as infile:
                for idx, (read_id, sequence, plus, quality) in enumerate(grouper(infile, 4, fillvalue=None)):
                    if read_id is None:
                        break

                    match = pattern.search(read_id)
                    if match:
                        found_read_id = match.group(1)
                        if found_read_id in readIDs:
                            outfile.write(f"{read_id}")
                            outfile.write(f"{sequence}")
                            outfile.write(f"{plus}")
                            outfile.write(f"{quality}")
    else:
        return


_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=XB])')
def _cigar_query_aln_len(cigar: str) -> int:
    if not cigar or cigar == "*":
        return 0
    qlen = 0
    for n, op in _CIGAR_RE.findall(cigar):
        n = int(n)
        # query-consuming ops; drop "S" to be stricter
        if op in ("M", "I", "=", "X"):
            qlen += n
    return qlen


def evaluate_bam_headerless(bam_path, thresholds):
    # dump primary-ish reads (you already filtered when creating the BAM,
    # but doing -F again is cheap/safe)
    proc = subprocess.run(
        ["samtools", "view", "-F", "4", "-F", "256", "-F", "2048", bam_path],
        text=True,
        capture_output=True
    )
    if proc.returncode != 0:
        print(f"[ERROR] samtools view {bam_path}: {proc.stderr}")
        return []

    lines = proc.stdout.strip().splitlines()
    total = len(lines)
    if total == 0:
        return []

    nm_vals = []
    aln_vals = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 11:
            continue
        cigar = parts[5]
        aln_len = _cigar_query_aln_len(cigar)
        if aln_len <= 0:
            continue

        nm = None
        for field in parts[11:]:
            if field.startswith("NM:i:"):
                nm = int(field[5:])
                break
        if nm is None:
            continue

        nm_vals.append(nm)
        aln_vals.append(aln_len)

    if not nm_vals:
        return []

    df = pd.DataFrame({"NM": nm_vals, "aln_len": aln_vals})

    p = Path(bam_path)
    seq_id = p.parts[-2] if len(p.parts) >= 2 else p.stem

    results = []
    for thr in thresholds:
        passing = (df["NM"] / df["aln_len"] < thr).sum()
        ratio = passing / len(df) # fixed from total
        results.append({
            "bam_path": bam_path,
            "seqID": seq_id,
            "threshold": thr,
            "total_reads": total,
            "reads_with_nm": len(df),
            "passing_reads": passing,
            "passing_ratio": round(ratio, 5),
        })
    return results


def edit_dist_thresholding(editdist_out, prim_map):
    thresholds = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.5]
    rows = evaluate_bam_headerless(prim_map, thresholds)
    if rows:
        pd.DataFrame(rows).to_csv(editdist_out, index=False)
        print(f"Saved {editdist_out}")


# run the analysis for each seqID
def run_analysis(deduplicate, save_dir, _str_, seqID, reference_genome, fq, asm5, sample, spp_name_clean):
    print(f"Running analysis for {spp_name_clean}, asm5: {asm5}")
    spp_name_clean = spp_name_clean.replace(" ", "_")
    if asm5:
        folder_out = f"{save_dir}/asm5/"
        sample_out = f"{folder_out}/{sample}/"
        seqID_asm5_out = f"{sample_out}/{seqID}/"
        sam_aligned = f"{seqID_asm5_out}/{spp_name_clean}-{_str_}-asm5-alignment.sam"
    else:
        folder_out = f"{save_dir}/map-ont/"
        sample_out = f"{folder_out}/{sample}/"
        seqID_asm5_out = f"{sample_out}/{seqID}/"
        sam_aligned = f"{seqID_asm5_out}/{spp_name_clean}-{_str_}-map-ont-alignment.sam"

    os.makedirs(folder_out, exist_ok=True)
    os.makedirs(sample_out, exist_ok=True)
    os.makedirs(seqID_asm5_out, exist_ok=True)

    bam_aligned = sam_aligned.replace(".sam", ".bam")
    sorted_path = bam_aligned.replace(".bam", "-sorted.bam")
    mpileup_file= sam_aligned.replace(".sam", "-mpileup.tsv")
    flag_file   = sam_aligned.replace(".sam", "-full-flagstat.txt")
    prim_map    = sorted_path.replace("-sorted.bam", "-primary-map.bam")
    txt_map     = prim_map.   replace(".bam", "-coverage.txt")
    editdist_out= sam_aligned.replace(".sam", "-editdist.tsv")

    if not os.path.isfile(sam_aligned):
    # if deduplicate or not os.path.isfile(sam_aligned):
        align_to_minimap2(fq, reference_genome, asm5, sam_aligned)
    if not os.path.isfile(bam_aligned):
    # if deduplicate or not os.path.isfile(bam_aligned):
        sam_to_bam(sam_aligned, bam_aligned)
    if not os.path.isfile(mpileup_file):
    # if deduplicate or not os.path.isfile(mpileup_file):
        sort_bam(bam_aligned, reference_genome, sorted_path, mpileup_file)
    if not os.path.isfile(flag_file):
    # if deduplicate or not os.path.isfile(flag_file):
        flagstat_run_and_save(sorted_path, flag_file)
    if not os.path.isfile(prim_map):
    # if deduplicate or not os.path.isfile(prim_map):
        primary_map(sorted_path, prim_map)
    if not os.path.isfile(txt_map):
    # if deduplicate or not os.path.isfile(txt_map):
        coverage(prim_map, txt_map)
    # if not os.path.isfile(editdist_out):
    if deduplicate or not os.path.isfile(editdist_out):
        edit_dist_thresholding(editdist_out, prim_map)
  
    return 


def process_task(task_data):
    """
    Function to process a single seqID task.
    """
    vs, viral_subset_df, fq_name, fq_seqID_unique, sample, directory = task_data

    spp_name = viral_subset_df.loc[viral_subset_df["seqID"] == vs]["name"].unique()[0]
    spp_name_clean = re.sub(r'[^\w\s]', '', spp_name)[0:40]
    readIDs = viral_subset_df.loc[viral_subset_df["seqID"] == vs]["readID"].unique()

    generate_fq_seqID_files(fq_name, fq_seqID_unique, sample, directory, vs, spp_name_clean, readIDs)
    return


def process_coverage_task(task_data):
    """
    Function to process a single coverage task.
    """
    # print(f"Number of variables: {len(task_data)}")
    fq_name, db_folder, database, fq_seqID_unique, full_fastq_bool, vs, asm5, sample, directory, combined_dict, viral_seqIDs, save_dir = task_data
    print(f"Processing: {vs}")
    try:
        if vs in viral_seqIDs:

            if asm5:
                seqID_asm5_out = f"{save_dir}/asm5/{sample}/{vs}/"
                mapped_file = f"{seqID_asm5_out}/*-{database}*asm5*-coverage.txt"
                if full_fastq_bool:
                    mapped_file = f"{seqID_asm5_out}/full*-{database}*asm5*-coverage.txt"
            else:
                seqID_asm5_out = f"{save_dir}/map-ont/{sample}/{vs}/"
                mapped_file = f"{seqID_asm5_out}/*-{database}*map-ont*-coverage.txt"
                if full_fastq_bool:
                    mapped_file = f"{seqID_asm5_out}/full*-{database}*map-ont*-coverage.txt"
            glob_map = glob.glob(mapped_file)
            if len(glob_map) > 0:
                mapped_file = glob_map[0]

            # if os.path.exists(mapped_file):
            #     return

            reference_genome_filename = combined_dict[vs]
            print(reference_genome_filename, vs, sample, mapped_file)
            if full_fastq_bool:
                fq_seqID_file = fq_name
                spp_name_clean = f"full_fq_{vs}"
            else:
                fq_seqID_file = f"{directory}/analysis/sample_data/{sample}/{fq_seqID_unique}/{vs}*.fastq"
                glob_fq_seqID_file = glob.glob(fq_seqID_file)
                if len(glob_fq_seqID_file) > 0:
                    fq_seqID_file = glob_fq_seqID_file[0]
                    # print(f"Fastq file at: {fq_seqID_file}; file size: {os.stat(fq_seqID_file).st_size/1e3:.0f} KB")
                    spp_name_clean = "_".join(fq_seqID_file.split("/")[-1].split("_")[1:]).replace(".fastq", "")
                else:
                    return

            reference_genome = f"/mnt/usersData/NCBI_RVDB/virus_RVDBs/{db_folder}/{reference_genome_filename}"
            print(f"\nProcessing: {spp_name_clean} --- Reference genome: {reference_genome} --- file size: {os.stat(reference_genome).st_size/1e3:.0f} KB")

            print(f"Checking for {mapped_file} --- {os.path.exists(mapped_file)}")

            if os.path.exists(reference_genome) and os.path.exists(fq_seqID_file):
                run_analysis(False, save_dir, f"{database}", vs, reference_genome, fq_seqID_file, asm5, sample, spp_name_clean)
                return
    except:
        print(f"Error processing {vs}")


def finish_setup(fq_seqID_unique, full_fastq_bool, directory, viral_subset_df, viral_seqIDs, sample):
    fq_name  = f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"

    if not full_fastq_bool:
        fq_gen_tasks = [
            (vs, viral_subset_df, fq_name, fq_seqID_unique, sample, directory)
            for vs in viral_seqIDs
        ]
        print(f"Number of tasks: {len(fq_gen_tasks)}")
        with mp.Pool(processes=20) as pool:
            pool.map(process_task, fq_gen_tasks)



def run_minmap2_after_setup(db_folder, model_name, full_fastq_bool, directory, save_dir, asm5, combined_dict, viral_seqIDs, sample_names):
    cov_run_tasks = []  # Initialize an empty list for all tasks

    for sample in sample_names:
        local_viral_seqIDs = viral_seqIDs[sample]
        fq_name = f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"

        # Add tasks for the current sample to the overall list
        intermediate_cov_run_tasks = [
            (fq_name, db_folder, database, fq_seqID_unique, full_fastq_bool, vs, asm5, sample, directory, combined_dict, list(set(list(local_viral_seqIDs))), save_dir)
            for vs in local_viral_seqIDs
        ]
        cov_run_tasks.append(intermediate_cov_run_tasks)
    # flatten the list
    cov_run_tasks = [item for sublist in cov_run_tasks for item in sublist]

    print(f"Number of tasks: {len(cov_run_tasks)}")
    n = 0
    for task in cov_run_tasks:

        for t in task:
            n += 1
            if isinstance(t, dict):
                print(f"Skipping item {n}")
                continue
            elif isinstance(t, list):
                print(f"Item {n}: {t[0]}")
            else:
                print(f"Item {n} in each variable: {t} --- {type(t)}")
        break

    with mp.Pool(processes=20) as pool:  # Adjust `processes` as needed
        pool.map(process_coverage_task, cov_run_tasks)
    print(f"Completed all tasks.")
    sys.exit(0)


def find_mapped_files(save_dir, full_fastq_bool, model_name, asm5, sample_names):
    data_files = []
    for sample in sample_names:
        if asm5:
            seqID_asm5_out = f"{save_dir}/asm5/{sample}/**/"
            mapped_file = f"{seqID_asm5_out}/*-{model_name}*asm5*-coverage.txt"
            if full_fastq_bool:
                mapped_file = f"{seqID_asm5_out}/full*-{model_name}*asm5*-coverage.txt"
        else:
            seqID_asm5_out = f"{save_dir}/map-ont/{sample}/**/"
            mapped_file = f"{seqID_asm5_out}/*-{model_name}*map-ont*-coverage.txt"
            if full_fastq_bool:
                mapped_file = f"{seqID_asm5_out}/full*-{model_name}*map-ont*-coverage.txt"
        glob_mapped_files = glob.glob(mapped_file)
        for gmf in glob_mapped_files:
            if os.path.exists(gmf):
                if os.stat(gmf).st_size > 0:
                    data_files.append(gmf)
    print(mapped_file)
    return data_files


def clean_read_bases(read_bases):
    """
    Remove common mpileup annotations (like indel markers, read start '^' and mapping quality,
    and read end '$') from the read bases string.
    This is a simplified cleaner; adjust according to your specific mpileup format.
    """
    # Remove the beginning of a read marker and its mapping quality (e.g., "^F")
    import re
    read_bases = re.sub(r'\^.', '', read_bases)
    # Remove end-of-read marker "$"
    read_bases = read_bases.replace("$", "")
    # Remove indel annotations: + or - followed by digits and bases.
    read_bases = re.sub(r'[\+\-]\d+[ACGTNacgtn]+', '', read_bases)
    return read_bases


def generate_consensus_from_mpileup(mpileup_file):
    consensus = []
    with open(mpileup_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue  # Skip if not enough columns
            ref_base = fields[2]
            coverage = int(fields[3])
            if coverage == 0:
                consensus.append('N')
            else:
                read_bases = clean_read_bases(fields[4])
                # Count occurrences of each base (ignoring case)
                base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                for base in read_bases.upper():
                    if base in base_counts:
                        base_counts[base] += 1
                # Determine consensus base using majority vote
                if sum(base_counts.values()) > 0:
                    consensus_base = max(base_counts, key=base_counts.get)
                else:
                    consensus_base = ref_base  # fallback
                consensus.append(consensus_base)
    return ''.join(consensus)


def calculate_sequence_identity(seq1, seq2):
    """
    Calculate sequence identity using global alignment.
    Uses Bio.Align.PairwiseAligner instead of pairwise2 for improved performance.

    :param seq1: First sequence (e.g., consensus sequence)
    :param seq2: Second sequence (e.g., reference sequence)
    :return: Identity percentage
    """

    # Defensive check
    if not seq1 or not seq2:
        print(f"⚠️ Skipping identity calc — one or both sequences empty "
              f"(len1={len(seq1)}, len2={len(seq2)})")
        return np.nan  # or 0.0 if you prefer a numeric fallback

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1

    alignment_iter = aligner.align(seq1, seq2)
    try:
        best_alignment = next(iter(alignment_iter))
    except StopIteration:
        return float("nan")

    # alignments = aligner.align(seq1, seq2)

    # if not alignments:
    #     print("⚠️ No alignment returned — sequences may be incompatible.")
    #     return np.nan

    # best_alignment = alignments[0]  # Best scoring alignment

    aligned_seq1 = best_alignment.aligned[0]
    aligned_seq2 = best_alignment.aligned[1]

    # Convert aligned segments into full aligned sequences
    aligned_seq1_str = ''.join(seq1[start:end] for start, end in aligned_seq1)
    aligned_seq2_str = ''.join(seq2[start:end] for start, end in aligned_seq2)

    # Compute identity percentage
    matches = sum(1 for a, b in zip(aligned_seq1_str, aligned_seq2_str) if a == b)
    alignment_length = len(aligned_seq1_str)

    if alignment_length == 0:
        return 0.0

    identity_percentage = (matches / alignment_length) * 100
    return identity_percentage


# --- helper just for edit-distance ---
def get_editdist_metrics(editdist_tsv: str, target_threshold: float) -> dict:
    """
    Given a ...-primary-map-coverage.txt file, try to load the corresponding
    ...-alignment-editdist.tsv and return the row for the chosen threshold.

    Always returns a dict with the same keys, filled with None on failure.
    """
    # # infer path
    # backup_editdist_file = editdist_tsv.replace("-uviral25-2", ".0-uviral25-2")
    # # in case your naming is just "-editdist.tsv" (older runs), try that too
    # fallback_file = mapfile.replace("-primary-map-coverage.txt", "-editdist.tsv")

    base_result = {
        "ED_threshold": None,
        "ED_total_reads": None,
        "ED_reads_with_nm": None,
        "ED_passing_reads": None,
        "ED_passing_ratio": None,
    }
    
    print(editdist_tsv)
    # print(backup_editdist_file)

    try:
        if os.path.exists(editdist_tsv):
            print(f"\nLoading edit-distance data from {editdist_tsv}\n")
            eddf = pd.read_csv(editdist_tsv)
        # elif os.path.exists(backup_editdist_file):
        #     eddf = pd.read_csv(backup_editdist_file)
        else:
            # nothing to load
            return base_result

        # make sure required columns exist
        req = {"threshold", "total_reads", "reads_with_nm", "passing_reads", "passing_ratio"}
        if not req.issubset(eddf.columns):
            return base_result

        row = eddf.loc[eddf["threshold"] == target_threshold]
        if row.empty:
            # no exact match → take the loosest one (last row)
            row = eddf.iloc[[-1]]

        row = row.iloc[0]
        print(row)
        return {
            "ED_threshold": float(row["threshold"]),
            "ED_total_reads": int(row["total_reads"]),
            "ED_reads_with_nm": int(row["reads_with_nm"]),
            "ED_passing_reads": int(row["passing_reads"]),
            "ED_passing_ratio": float(row["passing_ratio"]),
        }
    except Exception as e:
        # fail gracefully
        print(f"[editdist] warning: could not load edit-distance data for: {e}")
        return base_result
# --- end helper ---


def needs_ed_rerun(intermediate_out_json: str, current_threshold: float) -> bool:
    """
    Return True if ED metrics are missing/null and analysis should be rerun.
    Return False if ED metrics are present and valid.
    """
    if not os.path.exists(intermediate_out_json):
        return True  # nothing exists → must run

    try:
        with open(intermediate_out_json, "r") as f:
            data = json.load(f)

        # unwrap one level if needed
        if len(data) == 1:
            data = next(iter(data.values()))

        ed_fields = [
            "ED_threshold",
            "ED_total_reads",
            "ED_reads_with_nm",
            "ED_passing_reads",
            "ED_passing_ratio",
        ]

        for field in ed_fields:
            val = data.get(field, None)

            # treat missing, null, or NaN as invalid
            if val is None:
                return True
            if isinstance(val, float) and np.isnan(val):
                return True

        stored_threshold = float(data["ED_threshold"])
        if abs(stored_threshold - float(current_threshold)) > 1e-9:
            return True

        return False

    except Exception as e:
        print(f"[WARN] Could not validate ED fields in {intermediate_out_json}: {e}")
        return True  # fail safe → rerun


def extract_mapping_data(_type_, edit_distance_threshold, data_files, directory, model_name, save_dir, asm5, column_names, viral_seqIDs_files, db_folder):

    ref_location = f"/mnt/usersData/NCBI_RVDB/virus_RVDBs/{db_folder}/"

    mapping_dict = {}
    missing_sIDs = []
    for mapfile in data_files:
        mpileup_file = mapfile.replace("-primary-map-coverage.txt", "-mpileup.tsv")
        print(mpileup_file)
        # primary_mapped_file = mapfile.replace("-primary-map-coverage.txt", "-primary-map.bam")
        # flag_file = mapfile.replace("-primary-map-coverage.txt", "-full-flagstat.txt")
        sub_df = pd.read_csv(mpileup_file, sep='\t', header=None, names=column_names)

        mapfile = mapfile.replace("//", "/")
        sID = mapfile.split("/")[-2]
        sample = mapfile.split("/")[-3]

        basefolder = "/".join(mapfile.split("/")[0:-1])
        print(model_name)
        if asm5:
            primary_mapped_file = f"{basefolder}/*{model_name}-asm5-alignment-primary-map.bam"
            if _type_ == "TaxID_":
                primary_mapped_file = f"{basefolder}/*taxID*.0-{model_name}-asm5-alignment-primary-map.bam"
        else:
            primary_mapped_file = f"{basefolder}/*{model_name}-map-ont-alignment-primary-map.bam"
            if _type_ == "TaxID_":
                primary_mapped_file = f"{basefolder}/*taxID*.0-{model_name}-map-ont-alignment-primary-map.bam"
           
        print(primary_mapped_file)                
        primary_mapped_file = glob.glob(primary_mapped_file)
        if len(primary_mapped_file) > 0:
            primary_mapped_file = primary_mapped_file[0]

            ref_length = 0
            covered_bases = 0
            if os.path.isfile(mapfile):
                if os.stat(mapfile).st_size > 0:
                    with open(mapfile, "r", encoding="utf-8") as mf:

                        if asm5:
                            intermediate_out_json = f"{save_dir}/asm5/{sample}/{sID}/{_type_}intermediate_mapping_stats.json"
                        else:
                            intermediate_out_json = f"{save_dir}/map-ont/{sample}/{sID}/{_type_}intermediate_mapping_stats.json"

                        if not needs_ed_rerun(intermediate_out_json, edit_distance_threshold):
                            print(f"[SKIP] ED complete for {sample} {sID}")
                            # load json and append to mapping_dict
                            with open(intermediate_out_json, "r") as f:
                                data = json.load(f)
                            mapping_dict[f"{sample}_{sID}"] = data[f"{sample}_{sID}"]
                            
                            continue
                        else:
                            print(f"[RERUN] ED missing → recomputing {sample} {sID}")

                        text = mf.read()
                        region_length = re.search(r'\((\d+\.\d+\S*|\d+\S*)\)', text)

                        covered_bases_regex = r"Covered bases:\s+(\d+(\.\d+)?(Kbp|bp))"
                        covered_bases = re.search(covered_bases_regex, text, re.MULTILINE)

                        ref_length = region_length.group(1) if region_length else None
                        covered_bases = covered_bases.group(1) if covered_bases else None

                        ref_length = convert_to_base_pairs(ref_length)
                        covered_bases = convert_to_base_pairs(covered_bases)

                        region_length_precise = ref_length
                        # covered_bases_precise = covered_bases
                        total_cov_bases = covered_bases
                        # copy the primary mapped bam using subprocess to output file
                        prim_map_copy = primary_mapped_file.replace(".bam", "_copy.bam")
                        shutil.copyfile(primary_mapped_file, prim_map_copy)

                        if os.path.isfile(primary_mapped_file):
                            if os.stat(primary_mapped_file).st_size > 0:
                                print(f"\n{os.path.isfile(primary_mapped_file)}, {primary_mapped_file}, {os.stat(primary_mapped_file).st_size}")
                                region_length_precise = get_total_reference_length(primary_mapped_file)
                                covered_bases_precise, total_cov_bases = get_covered_bases(primary_mapped_file)

                            if int(ref_length) % 100 == 0:
                                if region_length_precise != 0:
                                    ref_length = region_length_precise

                            if int(covered_bases) < 1000:
                                if total_cov_bases != 0:
                                    covered_bases = total_cov_bases

                            elif int(covered_bases) % 100 == 0:
                                if total_cov_bases != 0:
                                    covered_bases = total_cov_bases

                        if os.stat(primary_mapped_file).st_size > 0:
                            reference_length = get_reference_length(prim_map_copy)
                            total_covered_bases, mean_coverage = calculate_coverage_and_mean(prim_map_copy, reference_length)  # Get total and mean coverage

                            stats = calculate_bam_stats(prim_map_copy, sID, sample)
                            stats['Reference Length'] = reference_length
                            stats['Covered Bases'] = total_covered_bases
                            stats['Mean Coverage'] = mean_coverage

                            longest_read_data = get_longest_read(prim_map_copy)
                            longest_read_id = get_longest_read_id(prim_map_copy)
                            if longest_read_data:
                                stats['Longest Read'] = longest_read_data.split()[0]
                                stats['Longest Read ID'] = longest_read_id.split()[0]
                                fastq_output_path = f"{prim_map_copy.rsplit('.', 1)[0]}_longest_read.fastq"
                                print(fastq_output_path)
                                fq_path = f"{directory}/analysis/sample_data/{sample}/trimmed/t*q"
                                fq_path = glob.glob(fq_path)
                                if len(fq_path) > 0:
                                    fq_path = fq_path[0]
                                    extract_read_from_fastq(fq_path, longest_read_id, fastq_output_path)
                            else:
                                stats['Longest Read'] = "None"
                                stats['Longest Read ID'] = "None"
                            print(stats)
                        else:
                            stats = {
                                'Reference Length': "None",
                                'Covered Bases': "None",
                                'Mean Coverage': "None",
                                'Longest Read': "None",
                                'Longest Read ID': "None"
                            }

                        # Average coverage
                        average_coverage = sub_df['coverage'].mean()

                        # Coverage uniformity
                        coverage_std_dev = sub_df['coverage'].std()
                        if average_coverage and not np.isnan(average_coverage):
                            coverage_cv = coverage_std_dev / average_coverage
                        else:
                            coverage_cv = float('nan')  # Assign NaN if invalid denominator

                        # Coverage depth
                        depth_1X = (sub_df['coverage'] >= 1).mean()
                        depth_5X = (sub_df['coverage'] >= 5).mean()
                        depth_10X = (sub_df['coverage'] >= 10).mean()
                        depth_20X = (sub_df['coverage'] >= 20).mean()

                        # Coverage gaps
                        coverage_gaps = (sub_df['coverage'] == 0).sum()

                        # Calculate the GC content coverage
                        sub_df['Is_GC'] = sub_df['ref_base'].isin(['G', 'C'])
                        gc_coverage = sub_df.loc[sub_df['Is_GC'], 'coverage'].mean()
                        # GC Coverage Bias
                        if average_coverage and not np.isnan(average_coverage):
                            gc_coverage_bias = gc_coverage / average_coverage
                        else:
                            gc_coverage_bias = float('nan')  # Assign NaN if invalid denominator

                        # coverage length
                        coverage_length = sub_df[sub_df['coverage'] > 0].shape[0]
                        print(f"Coverage length: {coverage_length}")

                        # region mean targeted coverage
                        if ref_length and ref_length > 0:
                            region_mean_targeted_coverage = coverage_length / ref_length
                        else:
                            region_mean_targeted_coverage = float('nan')

                        # % sequence identity
                        consensus_sequence = generate_consensus_from_mpileup(mpileup_file)
                        # print("Consensus sequence:", consensus_sequence)

                        run_seqID_percentages = True
                        try:
                            if sID not in viral_seqIDs_files:
                                raise KeyError(f"{sID} not found in viral_seqIDs_files!")

                        except KeyError as e:
                            print(f"⚠️ WARNING: {e}")
                            run_seqID_percentages = False
                            missing_sIDs.append(sID)

                        reference_sequence = None
                        identity_percentage = "missing from seqID mapping file"
                        if run_seqID_percentages:
                            try:
                                reference_sequence_file = viral_seqIDs_files[sample][sID]
                            except:
                                reference_sequence_file = viral_seqIDs_files[sID]
                            ref_path = f"{ref_location}{reference_sequence_file}"
                            print(f"Reference sequence file: {ref_path}")

                            if os.path.isfile(ref_path):
                                record = SeqIO.read(ref_path, "fasta")
                                reference_sequence = str(record.seq)

                                # ✅ Check both sequences before attempting alignment
                                len_cons = len(consensus_sequence)
                                len_ref = len(reference_sequence)

                                if len_cons == 0 or len_ref == 0:
                                    print(f"⚠️ Skipping identity calculation for {sID}")
                                    print(f"    Consensus length: {len_cons}")
                                    print(f"    Reference length: {len_ref}")
                                    print(f"    Reference path: {ref_path}\n")
                                    identity_percentage = None  # or np.nan if you want to preserve numeric structure

                                identity_percentage = calculate_sequence_identity(consensus_sequence, reference_sequence)
                                print(f"Identity Percentage: {identity_percentage:.2f}%\n\n")

                        editdist_out= primary_mapped_file.replace("-primary-map.bam", "-editdist.tsv")
                        ed_metrics = get_editdist_metrics(editdist_out, target_threshold=edit_distance_threshold)

                        mapping_dict[f"{sample}_{sID}"] = {
                            'Sample': sample,
                            'SeqID': sID,
                            'Reference Length': ref_length,
                            'Covered Bases': covered_bases,
                            'Mean Coverage': average_coverage,
                            'Coverage Std Dev': coverage_std_dev,
                            'Coverage CV': coverage_cv,
                            'Depth 1X': depth_1X,
                            'Depth 5X': depth_5X,
                            'Depth 10X': depth_10X,
                            'Depth 20X': depth_20X,
                            'Coverage Gaps': coverage_gaps,
                            'GC Content Coverage': gc_coverage,
                            'GC Coverage Bias': gc_coverage_bias,
                            'Region Mean Targeted Coverage': region_mean_targeted_coverage,
                            'Identity Percentage': identity_percentage,
                            'Coverage Length': coverage_length,
                            **stats,
                            **ed_metrics,
                        }

                        # Replace with this
                        intermediate_serializable_dict = {
                            f"{sample}_{sID}": {k: convert_to_serializable(v) for k, v in mapping_dict[f"{sample}_{sID}"].items()}
                        }
                        with open(intermediate_out_json, 'w') as f:
                            json.dump(intermediate_serializable_dict, f, indent=4)
                        print(f"File size: {os.stat(intermediate_out_json).st_size/1e3:.0f} KB")

    # save dictionary to json
    model_name = model_name.replace("*-", "")
    if asm5:
        mapping_json_save = f"{save_dir}/asm5/{_type_}{model_name}_mapping_stats.json"
        # missing_sIDs_json = f"{save_dir}/asm5/missing_sIDs.json"
    else:
        mapping_json_save = f"{save_dir}/map-ont/{_type_}{model_name}_mapping_stats.json"
        # missing_sIDs_json = f"{save_dir}/map-ont/missing_sIDs.json"

    print(f"Saving mapping stats to {mapping_json_save}")
    serializable_mapping_dict = {key: {k: convert_to_serializable(v) for k, v in value.items()} for key, value in mapping_dict.items()}

    with open(mapping_json_save, 'w') as f:
        json.dump(serializable_mapping_dict, f, indent=4)

    # Create proper DataFrame - dig through the double nesting
    rows = []
    for outer_key, outer_dict in serializable_mapping_dict.items():
        # Check if there's an extra level of nesting
        if outer_key in outer_dict:
            # The key is repeated as a nested key - extract the inner dict
            inner_dict = outer_dict[outer_key]
            rows.append(inner_dict)
        else:
            # No extra nesting - use the outer dict directly
            rows.append(outer_dict)

    df = pd.DataFrame(rows)
    print(df)
    print("DataFrame columns:", df.columns.tolist())
    print("DataFrame shape:", df.shape)

    # Save the DataFrame to a CSV file
    if asm5:
        mapping_csv_save = f"{save_dir}/asm5/{_type_}{model_name}_mapping_stats.csv"
    else:
        mapping_csv_save = f"{save_dir}/map-ont/{_type_}{model_name}_mapping_stats.csv" # TaxID_uviral25-2_mapping_stats
    print(f"Saving mapping stats to {mapping_csv_save}")
    df.to_csv(mapping_csv_save, index=False)

    # Print a sample of the DataFrame
    print("DataFrame Preview:")
    print(df.head(), df.shape)


def add_labels(model_name, label_file, asm5, save_dir):
    # open the DataFrame from a CSV file
    if asm5:
        mapping_csv_save = f"{save_dir}/asm5/{model_name}_mapping_stats.csv"
        labelled_save_csv = f"{save_dir}/asm5/{model_name}_mapping_stats_labelled.csv"
    else:
        mapping_csv_save = f"{save_dir}/map-ont/{model_name}_mapping_stats.csv" # TaxID_uviral25-2_mapping_stats
        labelled_save_csv = f"{save_dir}/map-ont/{model_name}_mapping_stats_labelled.csv"
    df = pd.read_csv(mapping_csv_save)

    seq2tax_file = "/mnt/usersData/RVDB2/RVDB2_seqID.map"
    seq2tax = pd.read_csv(seq2tax_file, sep="\t", header=None, names=["SeqID", "TaxID"])
    df = df.merge(seq2tax, on="SeqID", how="left")

    label_df = pd.read_csv(label_file) # col names = ["code", "name"]
    # add a column ("name") to df by merging on "SeqID" and "code" columns
    df = pd.merge(df, label_df, left_on="SeqID", right_on="code", how="left")
    # save to csv
    print(f"Saving labelled mapping stats to {labelled_save_csv}")
    df.to_csv(labelled_save_csv, index=False)

    # 5) OPTIONAL filter if ED/ID/depth exist
    needed = ["ED_passing_ratio", "Identity Percentage", "Depth 1X"]
    if all(col in df.columns for col in needed):

        # ✳️ TUNABLE THRESHOLDS ✳️
        ED_MIN_RATIO = 0.50      # keep rows where ED_passing_ratio >= this
        ED_MIN_READS = 1         # optional: require at least this many ED_passing_reads

        RESCUE_MIN_IDENTITY = 80.0    # rescue clause: high identity
        RESCUE_MIN_DEPTH1X = 0.90     # and decent depth


        def _passes(row):
            # ED fields
            ed_ratio = 0.0 if pd.isna(row["ED_passing_ratio"]) else float(row["ED_passing_ratio"])
            ed_reads = 0
            if "ED_passing_reads" in row and not pd.isna(row["ED_passing_reads"]):
                ed_reads = int(row["ED_passing_reads"])

            # other fields
            ident = 0.0 if pd.isna(row["Identity Percentage"]) else float(row["Identity Percentage"])
            d1x = 0.0 if pd.isna(row["Depth 1X"]) else float(row["Depth 1X"])

            # rule 1: ED is strong
            if ed_ratio >= ED_MIN_RATIO and ed_reads >= ED_MIN_READS:
                return True

            # rule 2: rescue good mappings hurt by ED
            if ident >= RESCUE_MIN_IDENTITY and d1x >= RESCUE_MIN_DEPTH1X:
                return True

            return False

        filtered_df = df[df.apply(_passes, axis=1)].copy()

        if asm5:
            filtered_path = labelled_save_csv.replace(".csv", "_filtered.csv")
        else:
            filtered_path = (
                f"{save_dir}/map-ont/Third-Pass-FULL-untargeted-"
                f"{model_name}_mapping_stats_labelled_filtered.csv"
            )

        print(f"Saving filtered mapping stats to {filtered_path}")
        filtered_df.to_csv(filtered_path, index=False)
    else:
        print("ED/ID/Depth columns not found — skipping filtered output.")

    sys.exit(0)


def parse_arguments():
    """
    Parses command line arguments.

    Returns:
    argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Identify and Extract Reads Script.")
    parser.add_argument('-d', '--directory', type=str, required=True, help='Path to the fq directory')
    parser.add_argument('-o', '--save-dir', type=str, required=True, help='Path to the save directory')
    parser.add_argument('-s', '--sample-config', type=str, required=True, help='Path to the sample config file')
    parser.add_argument('-a', '--asm5', action='store_true', help='Use asm5 for minimap2')
    parser.add_argument('-f', '--fastq', action='store_true', help='Use full fastq')
    parser.add_argument('-x', '--extract-fastq', action='store_true', help='Extract fastq files')
    parser.add_argument('-r', '--run-alignment', action='store_true', help='Run alignment')
    parser.add_argument('-l', '--label-species', action='store_true', help='Label species')
    parser.add_argument('-db', '--database', type=str, required=True, help='Path to the database')
    parser.add_argument("--edit-distance", dest="edit_distance_threshold", type=float, default=0.15, help="Edit distance threshold for filtering alignments (default: 0.15).")
    return parser.parse_args()


def main(fq_seqID_unique, edit_distance_threshold, model_name, database, ce_name, db_folder, label_file, full_fastq_bool, label_species, extract_fastq, run_alignment, directory, save_dir, asm5, sample_names, seqid_deanonymise, json_save):

    if label_species:
        add_labels(model_name, label_file, asm5, save_dir)

    # open json file
    print(f"Opening {json_save}")
    with open(json_save, 'r') as f:
        combined_dict = json.load(f)

    print(f"Opening {seqid_deanonymise}")
    seqid_df = pd.read_csv(seqid_deanonymise)

    viral_taxIDs = [11768,61673,133704,11856,28285,10376,12814,11250,10891]
    # need to find ALL viral taxIDs.
    vertebrate_taxIDs = [9606,9544,9685,9823,10090,30586]
    viral_seqIDs = {}
    viral_seqIDs_files = {}

    for s in sample_names:
        fq_seqID = f"{directory}/analysis/sample_data/{s}/{fq_seqID_unique}/"
        viral_subset_df, _, viral_seqIDs_sample, _ = cov_analysis(db_folder, ce_name, database, s, save_dir, fq_seqID, seqid_df, combined_dict, viral_taxIDs, vertebrate_taxIDs)
        # subset comnined_dict on viral_seqIDs
        viral_targets_dict = {k: combined_dict[k] for k in viral_seqIDs_sample}
        # subset seqid_df on viral_seqIDs
        viral_targets_seqIDs = seqid_df.loc[seqid_df["seqID"].isin(viral_seqIDs_sample)]
        
        print(f"\nViral targets: {len(viral_targets_dict)}")
        print(f"First 10 viral targets: {list(viral_targets_dict.keys())[0:10]}")
        print(f"Viral seqIDs: {len(viral_targets_seqIDs)}")
        print(f"First 10 viral seqIDs: {viral_targets_seqIDs.head()}")
        
        if extract_fastq:
            finish_setup(fq_seqID_unique, full_fastq_bool, directory, viral_subset_df, viral_seqIDs_sample, s)
        
        viral_seqIDs_sample = list(set(list(viral_seqIDs_sample)))
        viral_seqIDs[s] = viral_seqIDs_sample
        viral_seqIDs_files[s] = viral_targets_dict

    # save viral_seqIDs to json
    viral_seqIDs_save = f"{save_dir}/{model_name}_viral_seqIDs.json"
    print(f"Saving viral seqIDs to {viral_seqIDs_save}")
    with open(viral_seqIDs_save, 'w') as f:
        json.dump(viral_seqIDs_files, f, indent=4)
    print(f"Saved viral seqIDs to {viral_seqIDs_save}")

    if extract_fastq:
        sys.exit(0)

    if run_alignment:
        run_minmap2_after_setup(db_folder, model_name, full_fastq_bool, directory, save_dir, asm5, combined_dict, viral_seqIDs, sample_names)

    # find mapped files
    column_names = ['chrom', 'pos', 'ref_base', 'coverage', 'pileup_bases', 'base_quals']
    data_files = find_mapped_files(save_dir, full_fastq_bool, model_name, asm5, sample_names)
    extract_mapping_data("", edit_distance_threshold, data_files, directory, model_name, save_dir, asm5, column_names, viral_seqIDs_files, db_folder)
    add_labels(model_name, label_file, asm5, save_dir)


if __name__ == "__main__":

    args = parse_arguments()
    directory = args.directory
    save_dir = args.save_dir
    sample_config = args.sample_config
    asm5 = args.asm5
    full_fastq_bool = args.fastq
    extract_fastq = args.extract_fastq
    run_alignment = args.run_alignment
    label_species = args.label_species
    database = args.database

    edit_distance_threshold = args.edit_distance_threshold
    
    samples = pd.read_csv(sample_config)
    samples["sample"] = samples["date"].astype("str")+"_"+\
        samples["NA"]+"_"+\
        samples["strain"]+"_"+\
        samples["concentration_CFU"]+"_"+\
        samples["batch"].astype("str")+"_"+\
        samples["duration_h"].astype("str")
    sample_names = samples["sample"].unique()

    databases = {
        "uviral25-2": {
            "seqid_deanonymise": "/mnt/usersData/NCBI_RVDB/virus_RVDBs/fixed_URVDB_seqID2.map",
            "json_save": "/mnt/usersData/NCBI_RVDB/virus_RVDBs/U-RVDB-clean-fix2-index.json",
            "label_file": "/mnt/usersData/NCBI_RVDB/virus_RVDBs/fixed_URVDB_seqID2.map",
            "db_folder": "URVDB_split_fa",
            "ce_name": "mini-u-viral2",
            "fq_seqID_unique": "fq_seqID_uv25_2",
            "model_name":"uviral25-2"
        },
        # "c-viral29": {
        #     "seqid_deanonymise": "/mnt/usersData/NCBI_RVDB/virus_RVDBs/fixed_CRVDB_V29_seqID.map",
        #     "json_save": "/mnt/usersData/NCBI_RVDB/virus_RVDBs/C-RVDB-V29-index.json",
        #     "label_file": "/mnt/usersData/NCBI_RVDB/virus_RVDBs/fixed_CRVDB_V29_seqID.map",
        #     "db_folder": "CRVDB_V29_split_fa",
        #     "ce_name": "c-viral29",
        #     "fq_seqID_unique": "fq_seqID_cv29",
        #     "model_name":"c-viral29"
        # },
    }

    found_db = databases[database]
    seqid_deanonymise = found_db["seqid_deanonymise"]
    json_save = found_db["json_save"]
    label_file = found_db["label_file"]
    db_folder = found_db["db_folder"]
    ce_name = found_db["ce_name"]
    fq_seqID_unique = found_db["fq_seqID_unique"]
    model_name = found_db["model_name"]

    main(fq_seqID_unique, edit_distance_threshold, model_name, database, ce_name, db_folder, label_file, full_fastq_bool, label_species, extract_fastq, run_alignment, directory, save_dir, asm5, sample_names, seqid_deanonymise, json_save)
