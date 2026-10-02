#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 11:16:11 2023

@author: mangi
"""

import subprocess
import os
from typing import List
import glob


def run_nanofilt(quality: str, input_file: str, \
                 output_directory: str, output_file: str) -> None:
    """
    Runs NanoFilt on the input file with the specified quality,
    and writes the output to the output directory.

    :param quality: Quality threshold for NanoFilt.
    :param input_file: Path to the input FASTQ file.
    :param output_directory: Path to the output directory.
    """

    # Run the NanoFilt command
    print(output_file)
    with open(output_file, 'w') as outfile:
        subprocess.run(['NanoFilt', '-q', str(quality), input_file], stdout=outfile)


def process_files(directory: str, q: str, adjust_qscore: str, samples: List) -> None:
    """
    Process all FASTQ files in the specified directory with NanoFilt.

    :param directory: Path to the directory containing FASTQ files.
    """

    for n, s in enumerate(samples):
        file = f"{directory}/analysis/sample_data/{s}/trimmed/trimmed_{s}.fastq"
        output_file =  f"{directory}/analysis/sample_data/{s}/trimmed/q{q}.fastq"
        if adjust_qscore is not None:
            file = f"{directory}/analysis/sample_data/{s}/trimmed/t-adjusted-{adjust_qscore}Q.fastq"
            output_file =  f"{directory}/analysis/sample_data/{s}/trimmed/t-adj-q{q}.fastq"
        
        files = glob.glob(file)
        if len(files) > 0:
            files = files[0]
            if not os.path.isfile(output_file):
                run_nanofilt(q, files, directory, output_file)
