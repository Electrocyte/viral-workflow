#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:42:59 2022

@author: mangi
"""

import glob
import re
from pathlib import Path
from itertools import zip_longest


def grouper(iterable, n, fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def find_IDs(processed_fq: str, n: int, truncated_fq: str) -> list:
    with open(processed_fq, "r") as a, open(truncated_fq, "w") as aa:
        print(f"Saving truncated fastq file to: {truncated_fq}")
        for readID, sequence, plus, quality in grouper(a, 4, None):
            truncated_sequence = sequence[:n]
            truncated_quality = quality[:n]
            aa.write(readID)
            aa.write(truncated_sequence)
            aa.write("\n")
            aa.write(plus)
            aa.write(truncated_quality)
            aa.write("\n")
        print(f"Save complete for {truncated_fq}.\n")


def main(trim_reads_files, no_host_trim_reads_files, n) -> (list, list):

    trunc_trim = []
    for trf in trim_reads_files:
        truncated_fq = trf.replace("/trimmed_2", f"/trim_truncated_{n}bp_2")
        print(truncated_fq)
        print(not Path(truncated_fq).is_file())
        if not Path(truncated_fq).is_file():
            find_IDs(trf, n, truncated_fq)
        trunc_trim.append(truncated_fq)

    trunc_no_host = []
    print(no_host_trim_reads_files)
    for nhtrf in no_host_trim_reads_files:
        nh_truncated_fq = nhtrf.replace("/no_host_2", f"/no_host_truncated_{n}bp_2")
        print(nh_truncated_fq)
        print(not Path(nh_truncated_fq).is_file())
        if not Path(nh_truncated_fq).is_file():  # if false skip
            find_IDs(nhtrf, n, nh_truncated_fq)
        trunc_no_host.append(nh_truncated_fq)

    return trunc_trim, trunc_no_host


if __name__ == "__main__":
    directory = "/mnt/usersData/Viral_CHO/"
    # directory = "E:/SequencingData/Viral_CHO/"
    bespoke_analyses = f"{directory}/analysis/bespoke/"
    multiple_fastq = "D:/GitHub/SMART-CAMP/configs/viral_DNA_all2.txt"
    multiple_fastq = "/home/james/SMART-CAMP/configs/viral_DNA_all2.txt"
    samples = [line.rstrip("\n").split(",") for line in open(f"{multiple_fastq}")]
    samples = [item for sublist in samples for item in sublist]
    init_reads = True
    n = 30

    trim_reads_files = []
    no_host_trim_reads_files = []
    for sample in samples:
        fq_file_type_trim = "trim"
        reads_dir = f"{directory}/analysis/sample_data/{sample}/trimmed"
        reads_file = f"{reads_dir}/{fq_file_type_trim}*q"
        glob_reads_files = glob.glob(reads_file)
        if len(glob_reads_files) > 0:
            trim_reads_files.append(glob_reads_files[0])

        fq_file_type_nh = "no_host"
        reads_dir = f"{directory}/analysis/sample_data/{sample}/trimmed"
        reads_file = f"{reads_dir}/{fq_file_type_nh}*q"
        glob_reads_files = glob.glob(reads_file)
        if len(glob_reads_files) > 0:
            no_host_trim_reads_files.append(glob_reads_files[0])
    main(trim_reads_files, no_host_trim_reads_files, n)


# time ./pipeline/pre_trim_fq.py
