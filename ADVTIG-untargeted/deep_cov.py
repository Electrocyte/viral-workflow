#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import argparse
import glob
import json
import re
import os
from itertools import zip_longest
from typing import Iterable, Optional, Tuple
from multiprocessing import Pool, cpu_count
import sys
from run_coverage import run_analysis, find_mapped_files, extract_mapping_data, add_labels


def build_taxon_columns(
    df: pd.DataFrame,
    min_hit_len: int = 0,
    max_cov_frac_med: float = 0.0,
    max_short_frac: float = 1.0,
) -> pd.DataFrame:
    """
    From a raw centrifuge troubleshooting DF with columns like:
        readID, seqID, taxID, hitLength, queryLength, ...
    produce one row per taxID with the extra flags/bins we want.

    Returns a per-taxID DataFrame.
    """

    # normalize coverage cutoff if user passed "10" meaning 0.10
    if max_cov_frac_med > 1:
        max_cov_frac_med = max_cov_frac_med / 100.0

    # normalize taxID to string
    def _norm_taxid(v):
        if pd.isna(v):
            return None
        s = str(v).strip()
        if s.lower() == "unclassified":
            return "0"
        try:
            f = float(s)
            return str(int(f)) if f.is_integer() else s
        except Exception:
            return s

    df = df.copy()
    df["taxID_norm"] = df["taxID"].apply(_norm_taxid)
    df["hitLength"] = pd.to_numeric(df["hitLength"], errors="coerce")
    df["queryLength"] = pd.to_numeric(df["queryLength"], errors="coerce")

    per_tax_rows = []

    for taxid, sub in df.groupby("taxID_norm", dropna=False):
        # ----- basic counts -----
        n_assign = len(sub)
        n_unique_reads = sub["readID"].nunique() if "readID" in sub.columns else n_assign

        # ----- long vs short -----
        h = sub["hitLength"].dropna()
        n_good_assign = int((h >= min_hit_len).sum())
        short_only_flag = (n_good_assign == 0 and len(h) > 0)

        # fraction short, relative to *this* min_hit_len
        frac_short = float((h < min_hit_len).sum()) / float(len(h)) if len(h) else 0.0

        # ----- per-read coverage -----
        # collapse to per read: for each read take max(hitLength) and its queryLength
        per_read = (
            sub.groupby("readID", as_index=False)
               .agg(queryLength=("queryLength", "first"),
                    maxHit=("hitLength", "max"))
        )
        per_read["cov_frac"] = (per_read["maxHit"] / per_read["queryLength"]).where(per_read["queryLength"] > 0)

        # median covered fraction for this taxon
        cov_frac_median = per_read["cov_frac"].median() if len(per_read) else np.nan

        # ----- flags/bins -----
        bin_len = "short_only" if short_only_flag else "has_long"

        # breadth flag: fail if median coverage fraction is below cutoff
        low_breadth_flag = bool(
            pd.isna(cov_frac_median) or (cov_frac_median < max_cov_frac_med)
        )
        breadth_bin = "low_breadth" if low_breadth_flag else "has_breadth"

        # optional short-fraction filter — don't drop here, just record if it exceeds
        exceeds_short_frac = frac_short > max_short_frac

        per_tax_rows.append({
            "taxID": taxid,
            "n_assignments": n_assign,
            "n_unique_reads": n_unique_reads,
            "frac_hitLength_lt_min": frac_short,     # same as frac_hitLength_lt100, but parametric
            "cov_frac_median": cov_frac_median,
            "short_only_flag": short_only_flag,
            "bin": bin_len,
            "low_breadth_flag": low_breadth_flag,
            "breadth_bin": breadth_bin,
            "exceeds_max_short_frac": exceeds_short_frac,
        })

    return pd.DataFrame(per_tax_rows)


def apply_stringency(per_tax: pd.DataFrame, params: dict) -> pd.DataFrame:
    mhl  = params["mhl"]
    mcfm = params["mcfm"]
    msf  = params["msf"]

    # we already encoded these into columns in build_taxon_columns
    df = per_tax.copy()

    # start with everything
    mask = pd.Series(True, index=df.index)

    # 1) drop taxa that had ONLY short hits
    #    (this is the main effect of raising mhl)
    mask &= ~df["short_only_flag"]

    # 2) drop taxa that failed breadth
    if mcfm > 0:
        mask &= ~df["low_breadth_flag"]

    # 3) drop taxa where too large a fraction of hits were short
    #    (we stored the test as a bool already)
    if msf < 1.0:
        mask &= ~df["exceeds_max_short_frac"]

    return df[mask]


def post_metagenomic_filter(test_params, ce_name, barrage_o_params, save_dir, sample_names):

    save_filter = f"{save_dir}/analysis/TaxID_SeqID_Sample.csv"

    if not os.path.exists(save_filter) or test_params:
        all_filtered = []
        for sample in sample_names:
            print(f"\n{sample}")
            cent_out = f"{save_dir}/analysis/sample_data/{sample}/centrifuge/"
            troubleshooting = f"{cent_out}{sample}_{ce_name}_centrifuge_troubleshooting_report.tsv"
            if not os.path.exists(troubleshooting):
                print(f"Error: {troubleshooting} does not exist")
                return
            ts_df = pd.read_csv(troubleshooting, sep="\t")
            ts_df["Sample"] = sample
            print(f"Troubleshooting df: {len(ts_df['taxID'].unique())} taxa")

            if test_params:
                for stringency_type, param_list in barrage_o_params.items():
                    per_tax = build_taxon_columns(
                        ts_df,
                        min_hit_len=param_list["mhl"],
                        max_cov_frac_med=param_list["mcfm"],
                        max_short_frac=param_list["msf"],
                    )
                    filtered = apply_stringency(per_tax, param_list)

                    print(f"Testing: {stringency_type}---mhl: {param_list['mhl']} / mcfm: {param_list['mcfm']} / msf:{param_list['msf']}. Filtered metagenome aligned taxa df: {len(filtered)}")
            else:
                per_tax = build_taxon_columns(
                    ts_df,
                    min_hit_len=barrage_o_params["Middle___"]["mhl"],
                    max_cov_frac_med=barrage_o_params["Middle___"]["mcfm"],
                    max_short_frac=barrage_o_params["Middle___"]["msf"],
                )
                filtered = apply_stringency(per_tax, barrage_o_params["Middle___"])     

                filtered = (
                    filtered[["taxID"]]               # from per-tax table
                    .assign(Sample=sample)               # we already know sample
                )

                ts_min = ts_df[["seqID", "taxID"]].drop_duplicates()
                filtered["taxID"] = filtered["taxID"].astype(str)
                ts_min["taxID"] = ts_min["taxID"].astype(str)

                filtered = filtered.merge(ts_min, on="taxID", how="left")

                # rename to match main mapping CSV
                filtered = filtered.rename(
                    columns={
                        "seqID": "SeqID",
                        "taxID": "TaxID",
                    }
                )
                print(f"Filtered metagenome aligned taxa df: {len(filtered['taxID'].unique())}")
                all_filtered.append(filtered)

        if all_filtered:
            filtered_keys = pd.concat(all_filtered, ignore_index=True).drop_duplicates()
        else:
            filtered_keys = pd.DataFrame(columns=["Sample", "SeqID", "TaxID"])
        print(filtered_keys)
        print(f"Saving to: {save_filter}")
        filtered_keys.to_csv(save_filter, index=False)
    
    else:
        print(f"Pre-loading target taxa from: {save_filter}")
        filtered_keys = pd.read_csv(save_filter)
    return filtered_keys


def grouper(iterable: Iterable, n: int, fillvalue: Optional[str] = None) -> Iterable[Tuple]:
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def deduplicate_fastq(fq_in: str, overwrite: bool = False) -> str:
    """
    Remove duplicate 4-line FASTQ entries.
    Creates <fq_in>.dedup.fastq and <fq_in>.dedup.stats.json.
    Returns the path to the deduped FASTQ.
    """
    fq_out = f"{fq_in}.dedup.fastq"
    stats_out = f"{fq_in}.dedup.stats.json"
    if os.path.exists(fq_out) and not overwrite:
        return fq_out                     # already done

    seen = set()
    n_total = n_unique = 0
    id_regex = re.compile(r'^@(\S+)')     # grab first token after '@'

    with open(fq_in, "r") as fin, open(fq_out, "w") as fout:
        for hdr, seq, plus, qual in grouper(fin, 4):
            if hdr is None:               # EOF safeguard
                break
            n_total += 1

            m = id_regex.match(hdr)
            read_id = m.group(1) if m else hdr.strip()

            if read_id not in seen:
                seen.add(read_id)
                fout.writelines([hdr, seq, plus, qual])
                n_unique += 1

    stats = dict(
        total_records=n_total,
        unique_records=n_unique,
        duplicated_records=n_total - n_unique,
    )
    with open(stats_out, "w") as js:
        json.dump(stats, js, indent=2)

    print(f"[dedup] {fq_in}: removed {stats['duplicated_records']} "
          f"duplicates ({n_unique}/{n_total} unique)")
    return fq_out


def load_file(model_name, save_dir):
    labelled_save_csv = f"{save_dir}/map-ont/{model_name}_mapping_stats_labelled.csv"
    df = pd.read_csv(labelled_save_csv)
    return df


def get_top_6_combined(group):
    top_mapped = group.sort_values("Mapped reads", ascending=False).head(6)
    top_coverage = group.sort_values("Region Mean Targeted Coverage", ascending=False).head(6)
    combined = pd.concat([top_mapped, top_coverage])
    return combined.drop_duplicates()


def get_seqids_by_taxid(df, seqid_col="SeqID"):
    """
    Group by TaxID and return unique seqIDs per TaxID.

    Parameters:
        df (pd.DataFrame): Filtered DataFrame.
        seqid_col (str): Column name that holds the sequence IDs.

    Returns:
        dict: {TaxID: list of unique seqIDs}
    """
    taxid_to_seqids = (
        df.groupby("TaxID")[seqid_col]
        .apply(lambda x: sorted(set(x.dropna())))
        .to_dict()
    )
    return taxid_to_seqids


def process_job(job):
    try:
        deduplicate, save_dir, model_name, seqID, reference_genome, output_fq_path, asm5, sample, spp_name_clean = job
        # {spp_name_clean}-{model_name}-*-alignment.sam
        run_analysis(deduplicate, save_dir, model_name, seqID, reference_genome, output_fq_path, asm5, sample, spp_name_clean)
    except Exception as e:
        print(f"[ERROR] Failed on job {job}: {e}")


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
    parser.add_argument('-db', '--database', type=str, required=True, help='Path to the database')
    parser.add_argument('-x', '--extract-mapping-data-bool', action='store_true', help='Extract mapping data')
    parser.add_argument('-a', '--add-labels-bool', action='store_true', help='Add labels to the data')
    parser.add_argument('-dd','--deduplicate', action='store_true', default=False)
    parser.add_argument("-e", "--edit-distance-threshold", type=float, default=0.1, help="Edit distance threshold for filtering alignments.")
    return parser.parse_args()


def main(deduplicate, edit_distance_threshold, label_file, add_labels_bool, extract_mapping_data_bool, local_cache, json_save, model_name, save_dir, sample_names, fq_seqID_unique):
    df = load_file(model_name, save_dir)
    asm5=False
    test_params = False

    barrage_o_params = {
            "None_____":{"mhl":0,  "mcfm":0.0,  "msf": 1.0},
            "Middle___":{"mhl":50, "mcfm":0.05, "msf": 0.75}, # -mcfm "0.05" -msf "0.75" -mhl "50"
            "Stringent":{"mhl":100,"mcfm":0.15, "msf": 0.5}   # -mcfm "0.15" -msf "0.5" -mhl "100"
        }
    
    filtered_keys = post_metagenomic_filter(test_params, ce_name, barrage_o_params, save_dir, sample_names)
    
    first_pass_extra_filt = f"{save_dir}/analysis/MetaFilt-First-Pass-untargeted-uviral25-2_mapping_stats_labelled.csv"
    if not os.path.exists(first_pass_extra_filt):
        print(f"Number of hits from First PASS: {len(df)}")
        df["Sample"] = df["Sample"].astype(str)
        df["SeqID"]  = df["SeqID"].astype(str)
        df["TaxID"]  = df["TaxID"].astype(str)
        
        df = df.merge(
            filtered_keys[["Sample", "SeqID", "TaxID"]].drop_duplicates(),
            on=["Sample", "SeqID", "TaxID"],
            how="inner",
        )
        print(f"Number of hits from First PASS with metagenomic filter applied: {len(df)}")
        print(f"Saving to: {first_pass_extra_filt}")
        df.to_csv(first_pass_extra_filt, index=False)

    third_filtered_path = (
        f"{save_dir}/map-ont/Third-Pass-FULL-untargeted-"
        f"TaxID_{model_name}_mapping_stats_labelled_filtered.csv"
    ) # Third-Pass-FULL-untargeted-TaxID_uviral25-2_mapping_stats_labelled_filtered
    third_pass_extra_filt = f"{save_dir}/analysis/MetaFilt-Third-Pass-untargeted-uviral25-2_mapping_stats_labelled.csv"
    print(f"Checking for: {third_filtered_path}")
    if os.path.exists(third_filtered_path):    
        print(f"Checking for: {third_pass_extra_filt}")
        # if not os.path.exists(third_pass_extra_filt):    
        third_df = pd.read_csv(third_filtered_path)
        print(f"Number of hits from Third PASS: {len(third_df)}")
        third_df = third_df.merge(
            filtered_keys[["Sample", "SeqID", "TaxID"]].drop_duplicates(),
            on=["Sample", "SeqID", "TaxID"],
            how="inner",
        )
        print(f"Number of hits from Third PASS with metagenomic filter applied: {len(third_df)}")
        print(f"Saving to: {third_pass_extra_filt}")
        third_df.to_csv(third_pass_extra_filt, index=False)

    # Ensure columns are numeric and handle missing data
    df["Mapped reads"] = pd.to_numeric(df["Mapped reads"], errors="coerce").fillna(0)
    df["Region Mean Targeted Coverage"] = pd.to_numeric(df["Region Mean Targeted Coverage"], errors="coerce").fillna(0)

    # Group by Sample and TaxID
    df_group = df.groupby(["Sample", "TaxID"])
    top_hits = df_group.apply(get_top_6_combined).reset_index(drop=True)

    out_file = f"{save_dir}/map-ont/{model_name}_top6_mapped_and_coverage.csv"
    top_hits.to_csv(out_file, index=False)
    print(f"Saved top 6 (mapped + coverage) per (Sample, TaxID) to {out_file}")

    if not os.path.exists(local_cache):
        with open(json_save, 'r') as f:
            combined_dict = json.load(f)
        # generate the local cache
        seqid_list = top_hits["SeqID"].unique().tolist()
        subset_of_seqIDs = {
            k: v for k, v in combined_dict.items()
            if k in seqid_list
        }
        print(f"Saving local cache to {local_cache}")
        with open(local_cache, 'w') as f:
            json.dump(subset_of_seqIDs, f)
        sys.exit(0)

    else:
        print(f"Loading local cache from {local_cache}")
        with open(local_cache, 'r') as f:
            combined_dict = json.load(f)
    
    original_model_name = model_name
    
    if add_labels_bool:
        # model_name = f"TaxID_taxID-{model_name}"
        model_name = f"TaxID_{model_name}" # TaxID_uviral25-2_mapping_stats_labelled_filtered
        add_labels(model_name, label_file, asm5, save_dir)
        sys.exit(0)

    if extract_mapping_data_bool:
        column_names = ['chrom', 'pos', 'ref_base', 'coverage', 'pileup_bases', 'base_quals']
        model_name = f"taxID-*-{model_name}"
        data_files = find_mapped_files(save_dir, False, model_name, asm5, sample_names)
        print(f"Beginning to extract mapping data for {len(data_files)} files")
        extract_mapping_data("TaxID_", edit_distance_threshold, data_files, directory, original_model_name, save_dir, asm5, column_names, combined_dict, db_folder)
        sys.exit(0)

    jobs_to_run = []
    for sample in sample_names:
        fq_file_loc = f"{directory}/analysis/sample_data/{sample}/{fq_seqID_unique}/"
        print(f"\nProcessing sample: {sample}")
        print(f"Looking for fastq files in: {fq_file_loc}")
        taxid_to_seqids = get_seqids_by_taxid(df, seqid_col="SeqID")
        for taxid, seqids in taxid_to_seqids.items():
            # print(f"TaxID {taxid} has {len(seqids)} unique seqIDs")
            found_files = []
            for seqid in seqids:
                seqid_file = f"{fq_file_loc}/{seqid}*.fastq"
                found_file = glob.glob(seqid_file)
                if found_file:
                    found_files.append(found_file[0])
                else:
                    pass
            if len(found_files) == 0:
                # print(f"Warning: No fastq files found for seqID {seqid} in TaxID {taxid}")
                continue
            print(f"Number of found fastq files for TaxID {taxid}: {len(found_files)}")
            # concatenate the fastq files into 1 per taxid
            output_fq_path = f"{fq_file_loc}/taxid_{taxid}_combined.fastq"

            # Concatenate the FASTQ files
            if not os.path.exists(output_fq_path):
                with open(output_fq_path, "w") as outfile:
                    for fq_file in found_files:
                        with open(fq_file, "r") as infile:
                            outfile.write(infile.read())

            if os.path.exists(output_fq_path):
                if deduplicate:
                    output_fq_path = deduplicate_fastq(output_fq_path)

            # print(f"Written combined FASTQ for TaxID {taxid} to {output_fq_path}")

            sample_hits = top_hits.loc[top_hits["Sample"] == sample]
            taxid_sample_hits = sample_hits.loc[sample_hits["TaxID"] == taxid]
            # print(f"\n{sample} TaxID {taxid} hits:")
            # print(taxid_sample_hits[["SeqID", "Mapped reads", "Region Mean Targeted Coverage", "Depth 1X", "Identity Percentage", "Reference Length","name"]])

            for seqID in taxid_sample_hits["SeqID"]:
                if seqID not in combined_dict:
                    print(f"[warning] SeqID {seqID} missing from seqID mapping file — skipping.")
                    continue  # skip this entry safely   
                reference_genome_fa = combined_dict[seqID]
                reference_genome = f"/mnt/usersData/NCBI_RVDB/virus_RVDBs/{db_folder}/{reference_genome_fa}"
                # print(reference_genome_fa)
                name = reference_genome_fa.split("_")[1:]
                name = "_".join(name).split(".")[0]
                name = name.replace(" ", "_")
                spp_name_clean = f"{name}-taxID-{taxid}"
                jobs_to_run.append((deduplicate, save_dir, model_name, seqID, reference_genome, output_fq_path, asm5, sample, spp_name_clean))

    # multiprocess the jobs
    print(f"\nLaunching {len(jobs_to_run)} jobs using multiprocessing...")
    print(f"Using {cpu_count()} CPU cores.")
    cores_to_use = max(1, int(cpu_count() * 0.25))
    print(f"Using {cores_to_use} CPU cores.")

    with Pool(processes=cores_to_use) as pool:
        pool.map(process_job, jobs_to_run)

    print("All jobs completed.")


if __name__ == "__main__":

    args = parse_arguments()
    directory = args.directory
    save_dir = args.save_dir
    sample_config = args.sample_config
    extract_mapping_data_bool = args.extract_mapping_data_bool
    add_labels_bool = args.add_labels_bool
    deduplicate = args.deduplicate

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
            "model_name":"uviral25-2",
            "local_cache": f"{save_dir}/map-ont/uviral25-2_local_cache.json"
        }
    }

    found_db = databases[database]
    seqid_deanonymise = found_db["seqid_deanonymise"]
    json_save = found_db["json_save"]
    label_file = found_db["label_file"]
    db_folder = found_db["db_folder"]
    ce_name = found_db["ce_name"]
    fq_seqID_unique = found_db["fq_seqID_unique"]
    model_name = found_db["model_name"]
    local_cache = found_db["local_cache"]

    main(deduplicate, edit_distance_threshold,label_file, add_labels_bool, extract_mapping_data_bool, local_cache, json_save, model_name, save_dir, sample_names, fq_seqID_unique)
