#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 17:43:40 2020

@author: Mangifera
"""
# blastn output interpreter

import pandas as pd
import glob
import os
import re
import warnings

warnings.filterwarnings("ignore")


def human_ref_handler(row):
    compiled = re.compile(r"\|(\w+.\w+)")
    return compiled.findall(row)[0]


# changed "\|(\w+)" to "\|(.\w+)" to capture:
# acc|GENBANK|S70572.1|endogenous
# acc|GENBANK|S70572.1|{endogenous    <- this item
# acc|GENBANK|V01180.1|v-mos-Mo []
# acc|GENBANK|EF179138.1|X Triticosecale sp. ZXZ
def RVDB_ref_handler(row):
    compiled = re.compile(r"(\w+)\|(\w+)\|([\w\.\d]+)\|(.*[\.\-\w]+)")
    return compiled.findall(row)[0][2]


def clean_qseqid(row):
    compile_clean = re.compile(r"(\w+\-\w+\-\w+\-\w+\-\w+$)")
    if len(row) < 40:
        return row
    if len(row) > 40:
        cleaned = compile_clean.findall(row)
        if len(cleaned) > 0:
            len_row = len(cleaned[0])
            remove = len_row - 36
            return cleaned[0][remove:]


# handle high speed blastn outputs
# https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# cols = query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
#   5.  mismatch    number of mismatches
#   6.  gapopen     number of gap openings
#   7.  qstart      start of alignment in query
#   8.  qend        end of alignment in query
#   9.  sstart      start of alignment in subject
#  10.  send        end of alignment in subject
def hs_blastn_import(directory: str) -> pd.DataFrame:
    df = pd.read_csv(
        f"{directory}",
        sep="\t",
        header=None,
        names=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatches",
            "gap_opens",
            "q_start",
            "q_end",
            "s_start",
            "s_end",
            "evalue",
            "bitscore",
        ],
    )
    df = df.dropna()
    df.reset_index(drop=True)
    df.qseqid = df.qseqid.apply(clean_qseqid)
    row_indices = df.dropna()

    df = df.loc[row_indices.index].reset_index()
    df = df.drop(["index"], axis=1)

    if df["sseqid"].str.contains("ref")[0]:
        df["sseqid"] = df["sseqid"].apply(human_ref_handler)

    if df["sseqid"].str.contains("acc")[0]:
        df["sseqid"] = df["sseqid"].apply(RVDB_ref_handler)

    df = df.rename(columns={"qseqid": "read_id"})
    if not df.empty:
        return df


def split_name(df: pd.DataFrame) -> pd.DataFrame:
    try:
        df['name'] = df['name'].astype(str)
        df_split = pd.DataFrame(
            df.name.str.split(" ", 2).tolist(), columns=["genus", "species", "strain"]
        )
    except TypeError as e:
        print(f"Error: {e}")
        # Handle the error or return the dataframe in its current state or empty
        df_split = pd.DataFrame(columns=["genus", "species", "strain"])

    df_split["name"] = df_split[["genus", "species", "strain"]].apply(
        lambda row: " ".join(row.values.astype(str)), axis=1
    )
    return df_split


# collect species genus, species and strain
def import_deanonymised_nomenclature(directory):
    seqid_deanonymiser = pd.read_csv(
        directory, sep=",", header=None, names=["sseqid", "name"]
    )
    seqid_deanonymiser = seqid_deanonymiser.iloc[1:]

    seqid_deanonymised_names = split_name(seqid_deanonymiser)
    seqid_deanonymised_clean = seqid_deanonymiser.merge(
        seqid_deanonymised_names, on="name", how="inner"
    )

    seqid_deanonymised_clean = seqid_deanonymised_clean.drop_duplicates(
        subset=["sseqid"]
    )
    return seqid_deanonymised_clean


# import blastN out data, format it neatly
def import_blastN_out(
    sseqids: pd.DataFrame,
    blastn_dir: str,
    search_term: str,
    bN_name: str,
    no_host: str,
    sample: str,
) -> (list, list):
    summary_files = glob.glob(
        f"{blastn_dir}/*{sample}{no_host}_{bN_name}.{search_term}.tsv", recursive=True
    )
    if len(no_host) > 0:
        if len(summary_files) == 0:
            summary_files = glob.glob(
                f"{blastn_dir}/*_{bN_name}.{search_term}.tsv", recursive=True
            )
            no_host = ""
    print(summary_files, f"{blastn_dir}/*{sample}{no_host}_{bN_name}.{search_term}.tsv")
    df_summaries = []
    id_name = []
    for summary_file in summary_files:
        try:
            if os.stat(summary_file).st_size > 0:
                summary_file = summary_file.replace("\\", "/")
                print(f"File ({summary_file}) contains data.")
                import_file = hs_blastn_import(summary_file)
                if import_file is not None:
                    seqid_deanonymised = import_file.merge(
                        sseqids, on="sseqid", how="inner"
                    )
                    df_summaries.append(seqid_deanonymised)
                    if "\\" in summary_file:
                        summary_file = summary_file.replace("\\", "/")
                    id_ = summary_file.split("/")[-1].split(".")
                    id_name.append(id_[0])
            else:
                print(f"Empty file ({summary_file}).")
        except OSError:
            print("No file.")
    return df_summaries, id_name


# extract summary information from run
def generate_summary_data(read_info: list, sample: str) -> pd.DataFrame:
    print(f"Current sample: {sample}")
    try:
        describe_df = read_info.groupby(["genus", "species", "sseqid"])[
            [
                "length",
                "pident",
                "bitscore",
                "mismatches",
                "gap_opens",
                "evalue",
                "mean_qscore_template",
            ]
        ].describe()
    except:
        describe_df = read_info.groupby(["genus", "species", "sseqid"])[
            ["length", "pident", "bitscore", "mismatches", "gap_opens", "evalue"]
        ].describe()

    describe_df.columns = ["_".join(a) for a in describe_df.columns.to_flat_index()]
    describe_df.reset_index(inplace=True)
    describe_df["sample"] = sample
    return describe_df


def main(
    blastn_dir: str,
    sseqids: pd.DataFrame,
    plots: bool,
    search_term: str,
    BLASTn_name: str,
    no_host: str,
    sample: str,
    directory: str,
) -> int:

    skip_describe = False
    n = 0
    print(
        f"\nRunning blastN interpreter.\nInput parameters: \nblastn_dir {blastn_dir}, \nplots: {plots}, \nsearch_term: {search_term}, \nblastN_name: {BLASTn_name}"
    )
    try:
        fastq_file = f"{directory}/{sample}/**/sequencing_summ*"
        fastq_seq_summ_glob = glob.glob(fastq_file)
        print(f"Searching for qscores here: {fastq_seq_summ_glob}")
        qscores = pd.DataFrame()
        if len(fastq_seq_summ_glob) > 0:
            fssg = fastq_seq_summ_glob[0]
            if os.path.isfile(fssg):
                fastq_summary_df = pd.read_csv(fssg, delimiter="\t")
                qscores = fastq_summary_df[["read_id", "mean_qscore_template"]]
                print(f"\nqscores for {sample}: \n{qscores}\n")

        # print(f"Current working directory: {blastn_dir}; search term: {search_term}, host reads removed: {no_host}")
        hs_classified_reads = 0
        if len(os.listdir(blastn_dir)) > 0:
            df_summaries, id_name = import_blastN_out(
                sseqids, blastn_dir, search_term, BLASTn_name, no_host, sample
            )
            if len(df_summaries) > 0:
                BLASTn_out = f"{blastn_dir}{BLASTn_name}"
                os.makedirs(BLASTn_out, exist_ok=True)
                describe_BLAST_file = f"{BLASTn_out}/describe_all_predictions.csv"
                # print(df_summaries)
                read_info = df_summaries[0]

                if len(qscores) > 0:
                    read_info = read_info.merge(qscores, on="read_id")
                    n += 1
                    print(f"Successfully added qscores, samples completed: {n}")

                # for old files, will require turning this off
                describe_df = pd.DataFrame()
                if skip_describe:
                    if not os.path.isfile(describe_BLAST_file):
                        print(f"DF size: {len(read_info)}")
                        describe_df = generate_summary_data(read_info, sample)
                        print(f"\nID name: {id_name}")
                        describe_df["genus_species"] = describe_df[
                            ["genus", "species"]
                        ].agg(" ".join, axis=1)
                        df_concat_plot_reads = describe_df.nlargest(10, "length_count")
                        print(
                            f"Saving to: {BLASTn_out}/describe_all_predictions.csv; \n{BLASTn_out}/describe_top_ten.csv"
                        )
                        describe_df.to_csv(describe_BLAST_file, index=False)
                        df_concat_plot_reads.to_csv(
                            f"{BLASTn_out}/describe_top_ten.csv", index=False
                        )
                    else:
                        describe_df = pd.read_csv(describe_BLAST_file)

                    hs_classified_reads = describe_df.length_count.sum()

                else:
                    try:
                        print(f"DF size: {len(read_info)}")
                        describe_df = generate_summary_data(read_info, sample)
                        print(f"\nID name: {id_name}")
                        describe_df["genus_species"] = describe_df[
                            ["genus", "species"]
                        ].agg(" ".join, axis=1)
                        df_concat_plot_reads = describe_df.nlargest(10, "length_count")
                        print(
                            f"Saving to: {BLASTn_out}/describe_all_predictions.csv; \n{BLASTn_out}/describe_top_ten.csv"
                        )
                        describe_df.to_csv(describe_BLAST_file, index=False)
                        df_concat_plot_reads.to_csv(
                            f"{BLASTn_out}/describe_top_ten.csv", index=False
                        )
                        hs_classified_reads = describe_df.length_count.sum()
                    except:
                        print("Potential error: IndexError: list index out of range")

                return hs_classified_reads
        else:
            print("Empty file.")
    except OSError:
        print("No file.")


if __name__ == "__main__":
    # blastn_dir = "E:/SequencingData/Viral_CHO/analysis/sample_data/20210518_ssDNA_MVM10-1CHOK1Qiagen_500000000CFU_7_36/blastN/"
    # seqid_dir  = "D:/GitHub/SMART-CAMP/C_RVDB_seqids.csv"
    # sample = "20210518_ssDNA_MVM10-1CHOK1Qiagen_500000000CFU_7_36"
    # no_host = ""

    # samples_input = "/home/james/SMART-CAMP/configs/viral_DNA_all.txt"
    # samples_input = "D:/GitHub/SMART-CAMP/configs/viral_DNA_all.txt"
    samples_input = "/home/james/SMART-CAMP/configs/ALL_DNA_BAC_FUN.txt"
    # samples_input = "D:/GitHub/SMART-CAMP/configs/ALL_DNA_BAC_FUN.txt"
    if samples_input != None:
        samples = [line.rstrip("\n").split(",") for line in open(f"{samples_input}")]
        samples = [item for sublist in samples for item in sublist]

    seqid_dir = "/home/james/SMART-CAMP/C_RVDB_seqids.csv"
    # seqid_dir  = "D:/GitHub/SMART-CAMP/C_RVDB_seqids.csv"
    # seqid_dir  = "/home/james/SMART-CAMP/all_seqids.csv"
    # seqid_dir  = "D:/GitHub/SMART-CAMP/all_seqids.csv"
    sseqids = import_deanonymised_nomenclature(seqid_dir)

    plots = False
    search_Term = "hs_bn"  # example
    blastn_name = "cviral"
    # blastn_name = "filter_bacteria"
    # blastn_name = "fungal_all"
    no_host = "_no_host"
    directory = "/mnt/usersData/DNA/"
    # directory = "D:/SequencingData/Harmonisation/DNA/"

    for sample in samples:
        # blastn_dir = f"/mnt/usersData/Viral_CHO/analysis/sample_data/{sample}/blastN/"
        # blastn_dir = f"E:/SequencingData/Viral_CHO/analysis/sample_data/{sample}/blastN/"
        blastn_dir = f"/mnt/usersData/DNA/analysis/sample_data/{sample}/blastN/"
        # blastn_dir = f"D:/SequencingData/Harmonisation/DNA/analysis/sample_data/{sample}/blastN/"
        main(
            blastn_dir,
            sseqids,
            plots,
            search_Term,
            blastn_name,
            no_host,
            sample,
            directory,
        )
