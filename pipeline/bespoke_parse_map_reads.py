#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 11:29:58 2020

@author: Mangifera
"""

from pathlib import Path
import re
import pandas as pd
import glob
import os

from itertools import zip_longest

# saw this in the python docs. looks like exactly what you need
def grouper(iterable, n, fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


no_reads_comp = re.compile(r"Number of reads: ([0-9]+)")
covered_bases_comp = re.compile(r"Coveredbases:([\d]*[.]*[\d]*)([Mbp|Kbp|bp])")
percent_covered_comp = re.compile(r"Percentcovered:([\d]*[.]*[\d]*)")
mean_coverage_comp = re.compile(r"Meancoverage:([\d]*[.]*[\w]*[+]*[\w].*[^x])")
ref_genome_len_comp = re.compile(r"\(([\d]*[.]*[\d]*)([Mbp|Kbp|bp])")
Mean_baseQ_comp = re.compile(r"MeanbaseQ:([0-9]+[\.]\d+|\d+)")
Mean_mapQ_comp = re.compile(r"MeanmapQ:([0-9]+[\.]\d+|\d+|(\-\w+))")


def ID_re(
    s,
    date,
    NA,
    strain,
    CFU,
    batch,
    duration,
    sample_name,
    no_reads,
    line,
    covered_bases,
    percent_covered,
    mean_coverage,
    mean_baseQ,
    mean_mapQ,
):
    no_reads_results = no_reads_comp.findall(no_reads)[0]
    print(
        sample_name,
        no_reads,
        line,
        covered_bases,
        percent_covered,
        mean_coverage,
        mean_baseQ,
        mean_mapQ,
    )
    cb_results = covered_bases_comp.findall(covered_bases.replace(" ", ""))
    covered_bases_mapped = float(cb_results[0][0])
    if cb_results[0][1] == "K":
        covered_bases_mapped *= 1000
    if cb_results[0][1] == "M":
        covered_bases_mapped *= 1000000

    pc_results = percent_covered_comp.findall(percent_covered.replace(" ", ""))[0]
    mc_results = mean_coverage_comp.findall(mean_coverage.replace(" ", ""))[0]
    mc_results = mc_results.replace("x", "").replace("\n", "")

    rgl_results = ref_genome_len_comp.findall(sample_name.replace(" ", ""))[0]
    rgl_results_mapped = float(rgl_results[0])
    if rgl_results[1] == "K":
        rgl_results_mapped *= 1000
    if rgl_results[1] == "M":
        rgl_results_mapped *= 1000000

    mbQ_results = Mean_baseQ_comp.findall(mean_baseQ.replace(" ", ""))[0]
    mmQ_results = Mean_mapQ_comp.findall(mean_mapQ.replace(" ", ""))[0][0]

    results = [
        s,
        date,
        NA,
        strain,
        CFU,
        batch,
        duration,
        no_reads_results,
        covered_bases_mapped,
        pc_results,
        mc_results,
        rgl_results_mapped,
        mbQ_results,
        mmQ_results,
    ]
    print(results)
    return results


def find_IDs(directory: list, date, NA, strain, CFU, batch, duration) -> list:
    found_IDs = []
    # print(date, NA, strain, CFU, batch, duration)
    with open(directory, "r", encoding="utf8") as a:
        for (
            sample_name,
            no_reads,
            line,
            covered_bases,
            percent_covered,
            mean_coverage,
            mean_baseQ,
            mean_mapQ,
            line2,
            line3,
            line4,
            line5,
            line6,
        ) in grouper(a, 13, None):
            s = sample_name.split(" ")[0]
            print(s, sample_name, directory)
            if "human_input" in directory:
                if "NC_00" in s:
                    results = ID_re(
                        s,
                        date,
                        NA,
                        strain,
                        CFU,
                        batch,
                        duration,
                        sample_name,
                        no_reads,
                        line,
                        covered_bases,
                        percent_covered,
                        mean_coverage,
                        mean_baseQ,
                        mean_mapQ,
                    )

            else:
                results = ID_re(
                    s,
                    date,
                    NA,
                    strain,
                    CFU,
                    batch,
                    duration,
                    sample_name,
                    no_reads,
                    line,
                    covered_bases,
                    percent_covered,
                    mean_coverage,
                    mean_baseQ,
                    mean_mapQ,
                )

            found_IDs.append(results)
    return found_IDs


def extract_coverage_output(directory: str, sample: str) -> list:
    date, NA, strain, CFU, batch, duration = sample.split("_")
    if Path(directory).stat().st_size != 0:
        found_IDs = find_IDs(directory, date, NA, strain, CFU, batch, duration)
        return found_IDs
    else:
        return None


def process_coverage(
    directory: str,
    save_bespoke_analysis_directory_lib: str,
    sample: str,
    save_bespoke_analysis_directory: str,
    fq_file_type: str,
) -> (pd.DataFrame, pd.DataFrame):
    mapped_reads = f"{save_bespoke_analysis_directory_lib}*_mapped.txt"

    globbage = glob.glob(mapped_reads)
    df = pd.DataFrame()

    if len(globbage) > 0:
        # print(globbage)
        coverage_data = extract_coverage_output(globbage[0], sample)
        if coverage_data is not None:
            df = pd.DataFrame(coverage_data)
            df = df.rename(
                columns={
                    0: "species",
                    1: "date",
                    2: "NA",
                    3: "strain",
                    4: "concentration_CFU",
                    5: "batch",
                    6: "duration_h",
                    7: "no_reads_mapped",
                    8: "covered_bases(bp)",
                    9: "percent_covered",
                    10: "mean_coverage(X)",
                    11: "ref_genome_len(bp)",
                    12: "Mean baseQ",
                    13: "Mean mapQ",
                }
            )

    if df.empty:
        return None, None

    df = df.drop_duplicates(
        subset=[
            "species",
            "date",
            "NA",
            "strain",
            "concentration_CFU",
            "batch",
            "duration_h",
        ],
        keep="last",
    )
    df_mean = df.groupby(
        ["species", "date", "NA", "strain", "concentration_CFU", "batch", "duration_h"]
    ).mean()
    df_mean = df_mean.add_suffix("_mean")
    df_concat = pd.concat([df_mean], axis=1, join="outer")
    df.to_csv(
        f"{save_bespoke_analysis_directory}/{fq_file_type}_individual.csv", index=False
    )
    df_concat.to_csv(f"{save_bespoke_analysis_directory}/{fq_file_type}_aggregate.csv")

    return df, df_concat


def main(
    samples: list, directory: str, fq_file_type: str, species_name: str
) -> pd.DataFrame:
    dfs = []
    # print(samples, directory, fq_file_type, species_name)
    save_bespoke_analysis_directory_fq = ""
    for sample in samples:
        bespoke_analyses = f"{directory}/analysis/bespoke/"
        save_bespoke_analysis_directory = f"{bespoke_analyses}/{sample}/"
        save_bespoke_analysis_directory_fq = (
            f"{bespoke_analyses}/{sample}/{fq_file_type}/"
        )

        save_bespoke_analysis_directory_lib = (
            f"{save_bespoke_analysis_directory_fq}/{species_name}/"
        )
        df, df_concat = process_coverage(
            directory,
            save_bespoke_analysis_directory_lib,
            sample,
            save_bespoke_analysis_directory,
            fq_file_type,
        )
        if df is not None:
            dfs.append(df)

    print(f"Generating folder: {save_bespoke_analysis_directory_fq}")
    if len(dfs) > 0:
        os.makedirs(save_bespoke_analysis_directory, exist_ok=True)
        cat_df = pd.concat(dfs)
        print(
            f"Saving concatenated file to: {save_bespoke_analysis_directory}/cat_{fq_file_type}_individual.csv"
        )
        cat_df.to_csv(
            f"{save_bespoke_analysis_directory}/cat_{fq_file_type}_individual.csv",
            index=False,
        )
        return cat_df


if __name__ == "__main__":
    spp_location = "/home/james/SequencingData/Centrifuge_libraries/viral/library/refseq_viral/GCF_000850105.1_ViralProj14686_genomic_dustmasked.fna"
    directory = "/mnt/usersData/Viral_human/"
    fq_file_type_trim = "trim"
    fq_file_type_nh = "no_host"
    multiple_fastq = "/home/james/SMART-CAMP/configs/viral_DNA_all8.txt"
    # multiple_fastq = "D:/GitHub/SMART-CAMP/configs/viral_DNA_all9.txt"

    with open(spp_location, "r") as f:
        first_line = f.readline()
        species_name = "-".join(first_line.split(" ")[1:])
        species_name = species_name.replace(",", "")
        species_name = species_name.replace("$", "")
        species_name = species_name.replace("\n", "")

    samples = [line.rstrip("\n").split(",") for line in open(f"{multiple_fastq}")]
    samples = [item for sublist in samples for item in sublist]

    trim_cat_df = main(samples, directory, fq_file_type_trim, species_name)
    nh_cat_df = main(samples, directory, fq_file_type_nh, species_name)

    trim_cat_df["analysis"] = "trim"
    nh_cat_df["analysis"] = "no-host"

    nh_trim = [trim_cat_df, nh_cat_df]
    if len(nh_trim) > 1:
        cat_df_nh_trim = pd.concat(nh_trim)
        print(
            f"Saving trim-no-host concatenated file to: {directory}/analysis/bespoke/trim-no-host-{species_name}-aggregate.csv"
        )
        cat_df_nh_trim.to_csv(
            f"{directory}/analysis/bespoke/trim-no-host-{species_name}-aggregate.csv",
            index=False,
        )


# time ./pipeline/bespoke_parse_map_reads.py
