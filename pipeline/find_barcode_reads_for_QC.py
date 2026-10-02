#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 11:05:45 2022

@author: mangi
"""

import glob
import time
import re
import pandas as pd
from pathlib import Path

from itertools import zip_longest

# saw this in the python docs. looks like exactly what you need
def grouper(iterable, n, fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def find_IDs(processed_fq: str) -> list:
    found_IDs = []
    readID = re.compile(r"^@(\w.+)(\srunid)")
    with open(processed_fq, "r") as a:
        for line, nextline, nexterline, nextestline in grouper(a, 4, None):
            found_ID = readID.findall(line)
            if len(found_ID) > 0:
                found_IDs.append(found_ID[0][0])
    return found_IDs


def find_fastqs(ex_val: str, directory: str) -> list:
    print(f"{directory}/*{ex_val}*/**/fastq_pro*/c*.fastq")
    globs = glob.glob(f"{directory}/*{ex_val}*/**/fastq_pro*/c*.fastq")

    if len(globs) > 0:
        return globs

    elif len(globs) == 0:
        print(f"{directory}/*{ex_val}*/**/fastq_pro*/2*-T.fastq")
        globs = glob.glob(f"{directory}/*{ex_val}*/**/fastq_pass/2*-T.fastq")
        if len(globs) > 0:
            return globs
        else:
            return []
        

    else:
        return []


def find_fq_n_ids(ex_val, directory: str, ex_key: str):
    fastq_globs = find_fastqs(ex_val, directory)
    if len(fastq_globs) > 0:
        for fastq_glob in fastq_globs:
            sample = fastq_glob.replace("\\", "/").split("/")[-4:-3][0]
            print(f"Current sample: {sample}")
            output_directory = "/".join(fastq_glob.replace("\\", "/").split("/")[:-2])

            save_file = f"{output_directory}/sequencing_summary.txt"
            print(save_file)
            print(f"Sequencing summary location: {directory}/{ex_key}/**/seq*sum*t")
            if not Path(save_file).is_file():
                foundIDs = find_IDs(fastq_glob)
                ex_globs = glob.glob(f"{directory}/{ex_key}/**/seq*sum*t")
                if len(ex_globs) > 0:
                    seq_sum_df = pd.read_csv(ex_globs[0], delimiter="\t")
                    
                    print(seq_sum_df)
                    
                    sample_specific_reads = seq_sum_df.loc[
                        seq_sum_df["read_id"].isin(foundIDs)
                    ]
                    if not sample_specific_reads.empty:
                        
                        sample_specific_reads.to_csv(save_file, index=False, sep="\t")
    return fastq_globs


def run_extract_reads_for_qc(directory: str, experiments: dict):
    start = time.time()

    for ex_key, ex_val in experiments.items():
        print(ex_val)
        if isinstance(ex_val, str):
            find_fq_n_ids(ex_val, directory, ex_key)
        if isinstance(ex_val, list):
            for e_val in ex_val:
                find_fq_n_ids(e_val, directory, ex_key)

    end = time.time()

    print(f"Run time: {int(end) - int(start)} seconds")


if __name__ == "__main__":
    # directory = "E:/SequencingData/Viral_CHO/"
    directory = "/mnt/usersData/DNA/"
    experiments = {
        "S8B18": "_818_",
        "S8B19": "_819_",
        "S8B20": "_820_",
        "S8B21": "_821_",
    }
    # experiments = {"S8B15": "_815_", "S8B16": "_816_", "S8B17": "_817_"}
    run_extract_reads_for_qc(directory, experiments)

# ./pipeline/find_barcode_reads_for_QC.py
