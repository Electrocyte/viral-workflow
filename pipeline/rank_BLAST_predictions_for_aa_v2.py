#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 13:58:24 2021

@author: mangi
"""

import glob
import pandas as pd


# run from mp_comb_ce_BN_v3.py
def generate_BLAST_ranks(
    adventitious_agent: str, BLAST_file: str, directory: str, BLASTn_name: str
) -> pd.DataFrame:
    title_aa = adventitious_agent.title()
    lower_aa = adventitious_agent.lower()

    input_file = glob.glob(f"{BLAST_file}")
    print(input_file)
    df = pd.read_csv(input_file[0])
    df["sample"] = (
        df["date"].astype(str)
        + "_"
        + df["NA"]
        + "_"
        + df["strain"]
        + "_"
        + df["concentration_CFU"].astype(str)
        + "CFU_"
        + df["batch"].astype(str)
        + "_"
        + df["duration_h"].astype(str)
    )
    df["date"] = df["date"].astype(str)
    df["batch"] = df["batch"].astype(str)

    df_sorted = df.sort_values(
        [
            "date",
            "strain",
            "length_count",
        ],
        ascending=[True, True, False],
    ).reset_index(drop=True)

    ranks = []
    aa_vals = 0
    for sample in df_sorted["sample"].unique():
        df_sorted_intermediate = df_sorted.loc[df_sorted["sample"] == sample]
        date, na, strain, cfu, batch, time = sample.split("_")
        cfu = cfu.replace("CFU", "")
        isolate_df = df_sorted_intermediate.loc[
            (df_sorted_intermediate["date"] == date)
            & (df_sorted_intermediate["batch"] == batch)
            & (df_sorted_intermediate["strain"] == strain)
        ]
        isolate_df.reset_index(inplace=True, drop=True)
        isolate_df = isolate_df.drop_duplicates(
            subset=[
                "name",
                "date",
                "NA",
                "strain",
                "concentration_CFU",
                "batch",
                "duration_h",
            ],
            keep="first",
        )
        isolate_df.reset_index(inplace=True, drop=True)

        aa_vals = isolate_df[
            isolate_df["name"].str.contains(
                f"{title_aa}|{lower_aa}|{adventitious_agent}"
            )
        ].reset_index(drop=True)
        if len(aa_vals["length_count"]):
            aa_vals = aa_vals["length_count"][0]
        else:
            aa_vals = 0

        aa_local_indices = isolate_df[
            isolate_df["name"].str.contains(
                f"{title_aa}|{lower_aa}|{adventitious_agent}"
            )
        ].index
        total_indices = len(isolate_df.index)

        if len(aa_local_indices) > 0:
            rank = f"{str(int(aa_local_indices.values[0])+1)}/{str(total_indices)}"
            AA_rank = f"{str(int(aa_local_indices.values[0])+1)}"
        else:
            rank = "No target reads found"
            AA_rank = "No target reads found"
        rank_info = [
            date,
            na,
            strain,
            f"{cfu}CFU",
            batch,
            time,
            rank,
            AA_rank,
            total_indices,
            aa_vals,
        ]
        ranks.append(rank_info)

    df_ranks = pd.DataFrame(
        ranks,
        columns=[
            "date",
            "NA",
            "strain",
            "concentration_CFU",
            "batch",
            "duration_h",
            "rank",
            "AA",
            "no_ranks",
            "AA_reads",
        ],
    )
    df_ranks_sorted = df_ranks.sort_values(
        ["date", "rank"], ascending=[True, True]
    ).reset_index(drop=True)

    save_name = f"{directory}/describe_rank_analysis_for_{BLASTn_name}_{adventitious_agent}_from_BLAST_classification.csv"
    save_name = save_name.replace(" ", "_")
    print(save_name)
    df_ranks_sorted.to_csv(save_name, index=False)
    return df_ranks_sorted


if __name__ == "__main__":
    directory = "E:/SequencingData/Viral_CHO/analysis/"
    BLAST_name = "describe_cviral_no_host_all_agg.csv"
    filename = f"{directory}/{BLAST_name}"
    adventitious_agent = "Minute virus"
    BLASTn_name = "cviral"
    # directory = "D:/SequencingData/Harmonisation/DNA/analysis/"
    # BLAST_name = "hsblastn_fungal_all_no_host_all_agg.csv"
    # filename = f"{directory}/{BLAST_name}"
    # adventitious_agent = "Candida"

    df_ranks_sorted = generate_BLAST_ranks(
        adventitious_agent, filename, directory, BLASTn_name
    )
