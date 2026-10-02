#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 14:11:41 2022

@author: mangi
"""

import pandas as pd
import os


def get_nanostats_df(directory: str, samples: list) -> pd.DataFrame:
    dfs = []
    for sample in samples:
        nanostats = f"{directory}/analysis/sample_data/{sample}/summary-plots-log-transformed/NanoStats.txt"
        if os.path.isfile(nanostats):
            with open(nanostats, "r") as f:
                n = 0
                lines = {}
                for line in f:
                    n += 1
                    if n < 10:
                        line = line.replace("\n", "")
                        line = line.replace(" ", "")
                        split_lines = line.split(":")
                        a, b = split_lines
                        lines[a] = b
            nanostats_df = pd.DataFrame.from_dict(
                lines, columns=[f"{sample}"], orient="index"
            )
            nanostats_df_transposed = nanostats_df.T
            if not nanostats_df_transposed.empty:
                nanostats_df_transposed = nanostats_df_transposed.drop(
                    ["Generalsummary"], axis=1
                )
                dfs.append(nanostats_df_transposed)
    cat_nanostats = pd.concat(dfs)
    cat_nanostats.reset_index(inplace=True)
    print(
        f"Saving (cat_nanostats.py) to: {directory}/analysis/nanoplot_summary_data.csv"
    )
    cat_nanostats.to_csv(f"{directory}/analysis/nanoplot_summary_data.csv", index=False)
    print(cat_nanostats)
    return cat_nanostats


if __name__ == "__main__":
    # directory = 'E:/SequencingData/Viral_CHO/'
    directory = "/mnt/usersData/DNA/"
    # samples = ["20210420_ssDNA_MVM50-1CHOK1_2500000000CFU_3_42",
    #             "20210420_ssDNA_MVM100-1CHOK1_5000000000CFU_3_42",
    #             "20210420_ssDNA_MVM500-1CHOK1_25000000000CFU_3_42",
    #             "20210420_DNA_CHOK1_1000000CFU_3_42",
    #             "20210420_ssDNA_MVM_10000000000CFU_3_42",]
    samples_input = "/home/james/SMART-CAMP/configs/temp_aDNA_all5.txt"
    # samples_input = "D:/GitHub/SMART-CAMP/configs/CRAAM.txt"
    if samples_input != None:
        samples = [line.rstrip("\n").split(",") for line in open(f"{samples_input}")]
        samples = [item for sublist in samples for item in sublist]

    cat_nanostats = get_nanostats_df(directory, samples)

    # home_directory = "/mnt/e/SequencingData/CRAAM/"
    home_directory = "/mnt/usersData/DNA/"
    
    print(f"mkdir {home_directory}/Nanoplots/")
    for s in samples:
        print(f"mkdir {home_directory}/Nanoplots/{s}/")
        print(
            f"scp -r -P 1414 james@137.132.22.78:{directory}/analysis/sample_data/{s}/summary-plots-log-transformed/*  {home_directory}/Nanoplots/{s}"
        )
