#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd


def extract_rep(sample):
    # Example: '20230119_DNA_advtig-Rep1B-E50_10CFU_41_36' → 'Rep1B'
    parts = sample.split("-")
    for part in parts:
        if part.startswith("Rep"):
            return part
    return "ZZZ"  # fallback if not found


def taxid_read_filt(df, type_):
    taxid_read_max = df.groupby(["Sample", "TaxID"])["Total reads"].unique().reset_index()
    taxid_read_max["Total reads"] = taxid_read_max["Total reads"].apply(lambda x: x[0] if len(x) == 1 else x)

    genus_group = (
        df.groupby(["Sample", "TaxID"])["genus"]
        .unique()
        .reset_index()
    )
    taxid_read_max = taxid_read_max.merge(genus_group, on=["Sample", "TaxID"])
    seqID_group = (
        df.groupby(["Sample", "TaxID"])["SeqID"]
        .unique()
        .reset_index()
    )
    taxid_read_max = taxid_read_max.merge(seqID_group, on=["Sample", "TaxID"])

    seqID_group["SeqID_count"] = seqID_group["SeqID"].apply(len)
    seqID_group.drop(columns="SeqID", inplace=True)

    # Merge into main table
    taxid_read_max = taxid_read_max.merge(seqID_group, on=["Sample", "TaxID"])

    # order by Sample, taxid
    taxid_read_max["Rep"] = taxid_read_max["Sample"].apply(extract_rep)
    taxid_read_max = taxid_read_max.sort_values(["Rep", "TaxID"]).reset_index(drop=True)

    taxids = list(set(taxid_read_max["TaxID"].unique()))
    columns=['Rep', 'SeqID_count','TaxID', 'Total reads','Sample','genus', 'SeqID']
    taxid_read_max = taxid_read_max[columns]
    print(f"\n{taxid_read_max}\n")

    # Convert list columns to comma-separated strings
    taxid_read_max["SeqID"] = taxid_read_max["SeqID"].apply(lambda x: ", ".join(x) if isinstance(x, list) else x)
    taxid_read_max["genus"] = taxid_read_max["genus"].apply(lambda x: ", ".join(x) if isinstance(x, list) else x)

    save_csv = f"{path_dir}/{type_}_pass_taxid_summary.csv"
    print(f"\n{save_csv}")
    taxid_read_max.to_csv(save_csv, index=False)
    print(taxid_read_max.columns)
    return taxids


# Define path to your CSV file (update if needed)
path_dir = "/mnt/e/SequencingData/ADVTIG_primary/"
file_path1 = f"{path_dir}/FIRST_PASS_untargeted_uviral25_mapping_stats_labelled.csv"
file_path2 = f"{path_dir}/SECOND_PASS_untargeted_uviral25_mapping_stats_labelled_v2.csv"

# Read CSV with tab separator (since it's TSV-style)
df_taxid = pd.read_csv(file_path2)
df_taxid["genus"] = df_taxid["name"].str.split(" ",expand=True)[0] +" "+ df_taxid["name"].str.split(" ",expand=True)[1]

# Basic inspection
print("Columns:", df_taxid.columns.tolist())
print("Shape:", df_taxid.shape)

# Optional: display first few rows
print(df_taxid.head())

taxids = taxid_read_filt(df_taxid, "second")

df2 = pd.read_csv(file_path1)
df2["genus"] = df2["name"].str.split(" ",expand=True)[0] +" "+ df2["name"].str.split(" ",expand=True)[1]

df2_taxid = df2.loc[df2["TaxID"].isin(taxids)]
print(f"\n\n{df2_taxid}")

taxid_read_filt(df2_taxid, "first")
