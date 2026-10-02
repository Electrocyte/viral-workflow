#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 11:04:19 2023

@author: mangi
"""

import os
import pandas as pd
import glob
from itertools import zip_longest
import shutil
from typing import List


# saw this in the python docs. looks like exactly what you need
def grouper(iterable, n, fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def adjust_q_scores(line: str, adj_val: str) -> str:
    return ''.join([chr(ord(char) - int(adj_val)) for char in line.strip()]) + '\n'


def process_fastq(input_path: str, output_path: str, adj_val: str) -> None:
    output_dir = os.path.dirname(output_path)  # Get the directory part of output_path
    os.makedirs(output_dir, exist_ok=True)  # Create the directory if it does not exist

    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for readID, sequence, plus, quality in grouper(infile, 4, None):
            outfile.write(readID)
            outfile.write(sequence)
            outfile.write(plus)
            outfile.write(adjust_q_scores(quality, adj_val))


def process_samples(input_directory: str, samples: List, adj_val: str) -> None:

    for n, s in enumerate(samples):
        input_file_paths = glob.glob(f"{input_directory}/analysis/sample_data/{s}/trimmed/trimmed_{s}.fastq")
        print(f"Processing {len(input_file_paths)} FASTQ files for sample {s}...")
        if len(input_file_paths) > 0:
            input_file_path = input_file_paths[0]
            temp_output_file_path = f"{input_directory}/analysis/sample_data/{s}/trimmed/t-adjusted-{adj_val}Q.fastq"
            print(temp_output_file_path)
            if not os.path.isfile(temp_output_file_path):
                process_fastq(input_file_path, temp_output_file_path, adj_val)
