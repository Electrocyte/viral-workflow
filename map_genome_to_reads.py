#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 15:07:53 2021

@author: Mangifera
"""

import os
import re
import sys
import subprocess
import glob
import json
from pathlib import Path
import tempfile
import pandas as pd
# if there are issues, most likely due to modin.pandas
# import modin.pandas as pd
from multiprocessing import Pool
from functools import partial
import bespoke_parse_map_reads, pre_trim_fq
import argparse

encoding = "utf-8"

# assess reads by mapping them to appropriate genomes
# minimap2 -ax map-ont "/path/to/genome/chinese_hamster/" /path/to/sample.fastq > out_alignment.sam
def align_to_minimap2(
    save_minimap2: str, threads: str, full_fastq: str, output_dir: str, lib: str, RNA: bool, asm5: bool):
    print(f"Library used: {save_minimap2}, {lib}")
    if not Path(f"{output_dir}/{save_minimap2}.sam").is_file():
        print("\nRunning minimap2\n")
        print(f"RNA?: {RNA}")
        if RNA:
            # minimap2 -ax splice -uf -k14 ref.fa sample.fa
            # print("\n", "minimap2", "-ax", "splice", "-uf", "-k14", lib, full_fastq)
            full_fasta = full_fastq.replace("fastq", "fasta")
            print(f"\nfasta file: {full_fasta}")
            minimap2 = subprocess.run(
                ["minimap2", "-ax", "splice", "-uf", "-k14", lib, full_fasta], capture_output=True
            )
        elif asm5:
            minimap2 = subprocess.run(
                ["minimap2", "-ax", "asm5", lib, full_fastq], capture_output=True
            )
        else:
            minimap2 = subprocess.run(
                ["minimap2", "-ax", "map-ont", lib, full_fastq], capture_output=True
            )  # trimmed_ref_path = reference sequence ; demux_reads = reads to be aligned
        returncode = minimap2.returncode
        if returncode != 0:
            print("No alignment found")
        save_alignment = f"{output_dir}/{save_minimap2}.sam"
        print(f"\nSaving alignment to: {save_alignment}\n")
        with open(save_alignment, "w") as f:
            f.write(str(minimap2.stdout, encoding))


# convert SAM to BAM
# samtools view -b test_alignment.sam > test_alignment.bam
def sam_to_bam(sam_aligned: str, bam_aligned: str, threads: str, output_dir: str):
    print("Collecting aligned sam files to generate bam files...")
    globbed_sam = glob.glob(f"{output_dir}/*{sam_aligned}")
    for glob_sam in globbed_sam:
        bam_pre = glob_sam.split(sam_aligned)[:-1][0]
        S2B = subprocess.run(
            ["samtools", "view", "-b", glob_sam], capture_output=True
        )
        with open(f"{bam_pre}{bam_aligned}", "wb") as ff:
            ff.write(S2B.stdout)


# number of threads may cause error
# sort BAM
# samtools sort -o sorted.bam test_alignment.bam -@ 60
def sort_bam(bam_aligned: str, threads, output_dir: str, lib: str):
    print("Collecting aligned bam files for sorting...")
    globbed_bam = glob.glob(f"{output_dir}/*{bam_aligned}")
    for glob_bam in globbed_bam:
        bam_sort_pre = glob_bam.split(bam_aligned)[:-1][0]
        sorted_bam_post = "_sorted.bam"
        sorted_path = f"{bam_sort_pre}{sorted_bam_post}"
        subprocess.run(
            ["samtools", "sort", "-o", sorted_path, glob_bam, "-@", threads],
            capture_output=True,
        )

        mpileup_file = f"{output_dir}-mpileup.tsv"
        print(mpileup_file)
        mpile_out = subprocess.run(["samtools", "mpileup", "-f", lib, sorted_path], capture_output=True)
        mPIPE = str(mpile_out.stdout, encoding)
        with open(mpileup_file, "w") as g:
            g.write(mPIPE)
        print(f"{os.stat(mpileup_file).st_size/1e6:.0f} MB for {mpileup_file}")


# filter sorted bam files
# samtools view -b -F 0x04 -q 7 -o sorted_filtered.bam sorted.bam
def filter_bam(bam_sorted: str, threads: str, output_dir: str, primary: bool):
    print("Filtering sorted bam files...")
    globbed_filt_bam = glob.glob(f"{output_dir}/*{bam_sorted}")
    for glob_filt_bam in globbed_filt_bam:
        bam_sort_pre = glob_filt_bam.split(bam_sorted)[:-1][0]
        # filtered_bam_post = "_filtered.bam"
        filter_save = f"{bam_sort_pre}_sorted_filtered.bam"

        # Include: 4 (unmapped), 256 (secondary alignment), and 2048 (supplementary alignment)
        if not primary:
            subprocess.run(
                [
                    "samtools",
                    "view",
                    "-b",
                    "-F",
                    "0x04",
                    "-q",
                    "7",
                    "-o",
                    filter_save,
                    glob_filt_bam,
                ],
                capture_output=True,
            )

        # Exclude: 4 (unmapped), 256 (secondary alignment), and 2048 (supplementary alignment)
        else:
            subprocess.run(
                [
                    "samtools",
                    "view",
                    "-b",
                    "-F",
                    "0x904",
                    "-q",
                    "7",
                    "-o",
                    filter_save,
                    glob_filt_bam,
                ],
                capture_output=True,
            )

# number of threads may cause error
# index sorted & filtered bam file
# samtools index sorted.bam -@ 60
def index_bam(threads: str, output_dir: str) -> None:
    print("Indexing filtered and sorted bam files...")
    globbed_index = glob.glob(f"{output_dir}/*_sorted.bam")
    for glob_index_bam in globbed_index:
        subprocess.run(
            ["samtools", "index", glob_index_bam, "-@", threads], capture_output=True
        )


# mapping bam
# samtools view -b -F 4 sorted.bam > mapped.bam
def map_bam(bam_filt, bam_map, threads, output_dir) -> None:
    print("Mapping index bam files...")
    globbed_index = glob.glob(f"{output_dir}/*_sorted.bam")
    for glob_index in globbed_index:
        map_pre = glob_index.split("_sorted.bam")[:-1][0]
        print(os.stat(glob_index).st_size)
        print("samtools", "view", "-b", "-F", "4", glob_index)
        mapped = subprocess.run(
            ["samtools", "view", "-b", "-F", "4", glob_index], capture_output=True
        )
        with open(f"{map_pre}{bam_map}", "wb") as fff:
            fff.write(mapped.stdout)
        print(os.stat(f"{map_pre}{bam_map}").st_size)


# calculate genomic coverage
# generate output indication of successful mapping txt file generation
# samtools coverage -H -w 40 mapped.bam | tee mapped_reads.txt
def coverage(bam_map: str, txt_map: str, output_dir: str) -> None:
    print("Collecting mapped files to generate coverage analysis...")
    globbed_map = glob.glob(f"{output_dir}/*{bam_map}")
    for glob_map_bam in globbed_map:
        coverage_map = subprocess.run(
            ["samtools", "coverage", "-H", "-w", "40", glob_map_bam],
            capture_output=True,
        )
        PIPE = str(coverage_map.stdout, encoding)
        map_pre = glob_map_bam.split(bam_map)[:-1][0]
        with open(f"{map_pre}{txt_map}", "w") as ffff:
            ffff.write(PIPE)


def remove_extraneous_files(fastq_dir: str, extension) -> None:
    files_to_delete = glob.glob(f"{fastq_dir}/*{extension}", recursive=True)
    for file_path_del in files_to_delete:
        try:
            os.remove(file_path_del)
        except:
            print("Error while deleting file : ", file_path_del)


def extract_read_info_from_sam(sam_file):
    total_mapped_bases = 0  # Changed from total_bases
    mapped_lengths = []

    with open(sam_file, 'r') as f:
        for line in f:
            if not line.startswith("@"):  # If the line isn't a header
                parts = line.split("\t")
                seq = parts[9]

                # For mapped reads, the FLAG should not have the 4th bit set (i.e., 0x4).
                flag = int(parts[1])
                if not (flag & 0x4):
                    mapped_lengths.append(len(seq))
                    total_mapped_bases += len(seq)  # Only accumulate the length if the read is mapped

    average_mapped_length = sum(mapped_lengths) / len(mapped_lengths) if mapped_lengths else 0
    max_mapped_length = max(mapped_lengths) if mapped_lengths else 0

    return average_mapped_length, max_mapped_length, total_mapped_bases  # Changed from total_bases


def multiprocess_sample(
    directory: str,
    bespoke_analyses: str,
    threads_per_process: str,
    clean_up: bool,
    spp_location: str,
    species_name: str,
    fq_file_type: str,
    init_reads: bool,
    RNA: bool,
    asm5: bool,
    primary: bool,
    fq_location: str,
) -> None:

    sam_aligned = ".sam"
    bam_aligned = "_aligned.bam"
    bam_sorted = "_sorted.bam"
    bam_filt = "_filtered.bam"
    bam_map = "_mapped.bam"
    txt_map = "_mapped.txt"

    sample = fq_location.split("sample_data")[-1].split("trimmed")[0].replace("/", "")

    species_name = re.sub(r'[\\/:();*?"<>|]', '-', species_name)
    if len(species_name) > 30:
        species_name = species_name[0:10] +"-"+ species_name[-30:]
    print(f"\n{species_name}; {bespoke_analyses}/{sample}/\n")

    save_bespoke_analysis_directory = f"{bespoke_analyses}/{sample}/"
    save_bespoke_analysis_directory_fq = f"{bespoke_analyses}/{sample}/{fq_file_type}/"
    save_bespoke_analysis_directory_lib = (
        f"{save_bespoke_analysis_directory_fq}/{species_name}/"
    )
    print(save_bespoke_analysis_directory_fq)
    print(f"\nSave location: {save_bespoke_analysis_directory_lib}")

    os.makedirs(save_bespoke_analysis_directory, exist_ok=True)
    os.makedirs(save_bespoke_analysis_directory_fq, exist_ok=True)
    os.makedirs(save_bespoke_analysis_directory_lib, exist_ok=True)

    # Save the stats to a file
    cs = "_covstats.txt"
    stats_file_path = f"{save_bespoke_analysis_directory_lib}/{species_name}{cs}"

    if not os.path.isfile(
        f"{save_bespoke_analysis_directory_lib}/{species_name}{txt_map}"
    ):
        align_to_minimap2(
            species_name,
            threads_per_process,
            fq_location,
            save_bespoke_analysis_directory_lib,
            spp_location, RNA, asm5
        )

    print(stats_file_path, "\n\n")
    # if not os.path.isfile(stats_file_path):
    sam_path = f"{save_bespoke_analysis_directory_lib}/{species_name}{sam_aligned}"
    print(f"Sam path: {sam_path}")
    average_mapped_length, max_mapped_length, total_mapped_bases = extract_read_info_from_sam(sam_path)
    print(f"Average Mapped Length: {average_mapped_length}")
    print(f"Max Mapped Length: {max_mapped_length}")
    print(f"Total Mapped Bases: {total_mapped_bases}")  # Adjusted this print statement


    print(f"Coverage stats: {stats_file_path}")
    with open(stats_file_path, 'w') as f:
        f.write(f"Average Mapped Length: {average_mapped_length}\n")
        f.write(f"Max Mapped Length: {max_mapped_length}\n")
        f.write(f"Total Mapped Bases: {total_mapped_bases}\n")

    if not os.path.isfile(
        f"{save_bespoke_analysis_directory_lib}/{species_name}{txt_map}"
    ):
        threads_per_process = str(threads_per_process)
        sam_to_bam(
            sam_aligned,
            bam_aligned,
            threads_per_process,
            save_bespoke_analysis_directory_lib,
        )
        sort_bam(bam_aligned, threads_per_process, save_bespoke_analysis_directory_lib, spp_location)

        filter_bam(
            bam_sorted,
            threads_per_process,
            save_bespoke_analysis_directory_lib, primary)

        index_bam(threads_per_process, save_bespoke_analysis_directory_lib)
        map_bam(
            bam_filt, bam_map, threads_per_process, save_bespoke_analysis_directory_lib
        )
        coverage(bam_map, txt_map, save_bespoke_analysis_directory_lib)


######################################
# handle threads and processes
def threads_n_processes(processes: int) -> (int, str):
    processes = int(processes)
    threads_per_process = 1

    # control thread count
    if processes > 20:
        threads_per_process = int(processes / processes)

    if processes <= 20:
        if processes > 5:
            threads_per_process = int(processes / 5)

    if processes <= 5:
        if processes > 1:
            threads_per_process = int(processes * 5)

    processes_to_start = processes
    if processes_to_start > 20:
        processes_to_start = 20

    print(processes_to_start, threads_per_process)
    return processes_to_start, threads_per_process


######################################


# fastq_file_no is a text file listing sample names (should be same as actual sample folder)
def run_coverage(
    directory: str,
    processes: str,
    clean_up: bool,
    samples: list,
    bespoke_analyses: str,
    reads_files: list,
    spp_location: str,
    species_name: str,
    fq_file_type: str,
    init_reads: bool,
    RNA: bool,
    asm5: bool,
    primary: bool

) -> None:

    sam_aligned = ".sam"
    bam_aligned = "_aligned.bam"
    bam_sorted = "_sorted.bam"
    bam_filt = "_filtered.bam"
    bam_map = "_mapped.bam"

    processes_to_start, threads_per_process = threads_n_processes(processes)

    if "no_host" in fq_file_type and asm5:
        return

    func = partial(
        multiprocess_sample,
        directory,
        bespoke_analyses,
        threads_per_process,
        clean_up,
        spp_location,
        species_name,
        fq_file_type,
        init_reads,
        RNA,
        asm5,
        primary
    )

    with Pool(processes=processes_to_start) as p:
        p.map(func, reads_files)  # process data_inputs iterable with pool

    print("Multiprocessing and coverage analysis complete")
    for sample in samples:
        # save_bespoke_analysis_directory = f"{bespoke_analyses}/{sample}/"
        save_bespoke_analysis_directory_fq = (
            f"{bespoke_analyses}/{sample}/{fq_file_type}/"
        )
        save_bespoke_analysis_directory_lib = (
            f"{save_bespoke_analysis_directory_fq}/{species_name}/"
        )

        mapped_file = f"{save_bespoke_analysis_directory_lib}/*mapped.txt"
        file_glob = glob.glob(mapped_file)
        if len(file_glob) > 0:
            with open(file_glob[0], "r") as f_in, open(
                f"{save_bespoke_analysis_directory_lib}/short_list.txt", "w"
            ) as f_out:
                for line in f_in:
                    if ">" not in line:
                        seqid = line.split(" ")[0]
                        f_out.write(seqid)
                        # print(seqid)
                    if "Number of reads" in line:
                        # no_reads = line.split(" ")[-1]
                        f_out.write(line)

        if clean_up:
            remove_extraneous_files(save_bespoke_analysis_directory_lib, sam_aligned)
            remove_extraneous_files(save_bespoke_analysis_directory_lib, bam_aligned)
            remove_extraneous_files(save_bespoke_analysis_directory_lib, bam_sorted)
            remove_extraneous_files(save_bespoke_analysis_directory_lib, bam_filt)
            remove_extraneous_files(save_bespoke_analysis_directory_lib, bam_map)
            remove_extraneous_files(save_bespoke_analysis_directory_lib, ".bam.bai")

    # CHO_coverage_aggregator.centrifuge_to_coverage(data_output, top_dir, cov_name)


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def multiples(m, count):
    for i in range(count):
        print(i * m)


def get_reference_length_from_bam(bam_file_path: str) -> int:
    command = f"samtools idxstats {bam_file_path}"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.stdout:
        # Sum the lengths of all reference sequences
        total_length = sum(int(line.split('\t')[1]) for line in result.stdout.strip().split('\n') if line)
        return total_length
    return 0

def get_breadth_of_coverage_from_bam(bam_file_path: str) -> int:
    # Command to calculate breadth of coverage
    command = f"samtools depth {bam_file_path} | awk '($3>0) {{count++}} END {{print count}}'"
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.stdout.strip():
        # Convert to int if there's output
        return int(result.stdout.strip())
    else:
        # Return 0 if there's no output
        return 0




def main(
    samples: list,
    spp_location: str,
    directory: str,
    clean_up: bool,
    species_name: str,
    bespoke_analyses: str,
    processes: int,
    read_max_len: int,
    init_reads: bool,
    RNA: bool,
    asm5: bool,
    primary: bool,
    hum_fq: bool,
    other: str
) -> None:
    trim_reads_files = []
    no_host_trim_reads_files = []

    fq_file_type_trim = "unknown"
    for sample in samples:

        if hum_fq:
            fq_file_type_trim = "human"
        elif isinstance(other, str):
            fq_file_type_trim = other
        else:
            fq_file_type_trim = "trim"
        reads_dir = f"{directory}/analysis/sample_data/{sample}/trimmed"
        reads_file = f"{reads_dir}/{fq_file_type_trim}*q"
        print(f"Searching for reads here: {reads_file}")
        glob_reads_files = glob.glob(reads_file)
        print(f"Found target files: {len(glob_reads_files)}")
        if len(glob_reads_files) > 0:
            trim_reads_files.append(glob_reads_files[0])

        fq_file_type_nh = "no_host"
        reads_dir = f"{directory}/analysis/sample_data/{sample}/trimmed"
        reads_file = f"{reads_dir}/{fq_file_type_nh}*q"
        glob_reads_files = glob.glob(reads_file)
        if len(glob_reads_files) > 0:
            no_host_trim_reads_files.append(glob_reads_files[0])

    if init_reads:
        trunc_trim, trunc_no_host = pre_trim_fq.main(
            trim_reads_files, no_host_trim_reads_files, read_max_len
        )
        trim_reads_files = trunc_trim
        no_host_trim_reads_files = trunc_no_host
        fq_file_type_trim = "trim_truncate"
        fq_file_type_nh = "no_host_truncate"

    if isinstance(other, str) and asm5:
        fq_file_type_trim = f"asm5-{other}"
    elif not isinstance(other, str) and asm5:
        fq_file_type_trim = "asm5"
    elif not isinstance(other, str) and primary:
        fq_file_type_trim = "primary"

    species_name = re.sub(r'[\\/:();*?"<>|]', '-', species_name)
    if len(species_name) > 30:
        species_name = species_name[0:10] +"-"+ species_name[-30:]
    print(f"\nSpecies: {species_name}; {fq_file_type_trim} [Line 521]\n")

    run_coverage(
        directory,
        processes,
        clean_up,
        samples,
        bespoke_analyses,
        trim_reads_files,
        spp_location,
        species_name,
        fq_file_type_trim,
        init_reads,
        RNA,
        asm5,
        primary
    )
    trim_cat_df = bespoke_parse_map_reads.main(
        samples, directory, fq_file_type_trim, species_name
    )

    full_bases = []
    for sample in samples:
        bam_file_str = f"{directory}/analysis/bespoke/{sample}/trim/{species_name}/*_sorted_filtered.bam"
        bam_files = glob.glob(bam_file_str)
        for bam_file_path in bam_files:
            ref_length = get_reference_length_from_bam(bam_file_path)
            covered_bases = get_breadth_of_coverage_from_bam(bam_file_path)
            full_bases.append({
                "sample": sample,
                "species_name": species_name,
                "reference_length": ref_length,
                "covered_bases": covered_bases
            })
    df = pd.DataFrame(full_bases)
    save_full = f"{directory}/analysis/bespoke/full/"
    os.makedirs(save_full, exist_ok=True)

    n = 1
    json_base_path = f"{save_full}/full_bases"
    while os.path.exists(f"{json_base_path}_{n}.json"):
        n += 1
    json_path = f"{json_base_path}_{n}.json"
    csv_path = f"{save_full}/full_bases_{n}.csv"

    df.to_csv(csv_path, index=False)
    with open(json_path, 'w') as json_file:
        json.dump(full_bases, json_file, indent=4)
    print(f"Saved full bases data to {json_path}; {csv_path}")

    run_coverage(
        directory,
        processes,
        clean_up,
        samples,
        bespoke_analyses,
        no_host_trim_reads_files,
        spp_location,
        species_name,
        fq_file_type_nh,
        init_reads,
        RNA,
        asm5,
        primary
    )
    nh_cat_df = bespoke_parse_map_reads.main(
        samples, directory, fq_file_type_nh, species_name
    )


    species_name = species_name.replace("|", "")

    if "_" in fq_file_type_trim:
        fq_file_type_trim = fq_file_type_trim.replace("_", "-")

    print(f"\n{fq_file_type_trim} {trim_cat_df} [Line 594]")
    if trim_cat_df is None:
        sys.exit()
    trim_cat_df["analysis"] = fq_file_type_trim

    trim_cat_df.to_csv(f"{directory}/analysis/bespoke/trim-only-{species_name}-aggregate.csv",index=False)

    nh_trim = pd.DataFrame()
    if nh_cat_df is not None:
        nh_cat_df["analysis"] = fq_file_type_nh.replace("_", "-")
        nh_trim = [trim_cat_df, nh_cat_df]

    if len(nh_trim) > 1:
        cat_df_nh_trim = pd.concat(nh_trim)

        if not init_reads:
            print(
                f"Saving trim-no-host concatenated file to: {directory}/analysis/bespoke/trim-no-host-{species_name}-aggregate.csv"
            )
            cat_df_nh_trim.to_csv(
                f"{directory}/analysis/bespoke/trim-no-host-{species_name}-aggregate.csv",
                index=False,
            )

        if init_reads:
            print(
                f"Saving trim-no-host concatenated file to: {directory}/analysis/bespoke/trim-no-host-{species_name}-truncated-aggregate.csv"
            )
            cat_df_nh_trim.to_csv(
                f"{directory}/analysis/bespoke/trim-no-host-{species_name}-truncated-aggregate.csv",
                index=False,
            )


if __name__ == "__main__":
    # directory = "/mnt/usersData/Viral_CHO/"
    # multiple_fastq = "/home/james/SMART-CAMP/configs/viral_DNA_all2.txt"

    # directory = "/mnt/usersData/Viral_human/"
    # multiple_fastq = "/home/james/SMART-CAMP/configs/viral_DNA_all10.txt"

    parser = argparse.ArgumentParser(description='Running Time to Read')
    parser.add_argument('-d', '--directory', help='Directory with experiment.')
    parser.add_argument('-t', '--spp', help='Target species')
    parser.add_argument('-s', '--samples', help='Directory for sample names')
    parser.add_argument('-r', '--RNA', action='store_true', help='RNA is used', default=False)
    parser.add_argument('-a', '--asm5', action='store_true', help='Use asm5 instead of map-ont', default=False)
    parser.add_argument('-y', '--primary', action='store_true', help='Use only primarly alignment in map-ont', default=False)
    parser.add_argument('-f', '--hum_fq', action='store_true', help='Human reads only', default=False)
    parser.add_argument('-p', '--processes', help='Number of threads', default = 5)
    parser.add_argument('-o', '--other', help='Other files.')

    args = parser.parse_args()

    directory = args.directory
    input_type = args.spp
    sample_names = args.samples
    rna = args.RNA
    asm5 = args.asm5
    primary = args.primary
    hum_fq = args.hum_fq
    processes = args.processes
    other = args.other

    if "advtig-ebv-double-v2" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/EBV_V01555.5_doubled.fa"
    if "advtig-ebv-split-v2" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/EBV_V01555.6_splitted_v2.fa"
    if "advtig-pcv1-double-v2" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/PCV1_NC_001792.5_doubled_v2.fa"
    if "advtig-pcv1-split-v2" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/PCV1_NC_001792.6_splitted.fa"

    if "advtig-felv" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/FeLV.fa"
    if "advtig-adv5" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/ADV5.fa"
    if "advtig-felv2" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/FeLV_Kawakami-Theilen_strain.fa"
    if "advtig-pcv1" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/PCV1.fa"
    if "advtig-pcv1-double" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/fixed-PCV1-double.fa"
    if "advtig-pcv1-flip" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/NC_001792_Flipped.fa"
    if "advtig-reo" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REOvirus segments.fa"
    if "advtig-ebv" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/EBV-B95-8.fa"
    if "advtig-ebv-flip" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/fixed-EBV-B95-8-flipped.fa"
    if "advtig-rsv" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/RSV-A.fa"
    if "advtig-cat" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/Cat_GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fa"
    if "advtig-human" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/Human_GCF_000001405.40_GRCh38.p14_genomic.fa"

    if "advtig-primary-human" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/Human_hg38_p14_Primary_Assembly.fasta"

    if "advtig-tamarin" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/GCA_021498475.1_ASM2149847v1_genomic.fna"

    if "advtig-monkey" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/Monkey_GCF_003339765.1_Mmul_10_genomic.fa"
    if "advtig-pig" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/Pig_GCF_000003025.6_Sscrofa11.1_genomic.fa"
    if "advtig-reol1" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REO1_L1.fasta"
    if "advtig-reol2" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REO1_L2.fasta"
    if "advtig-reol3" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REO1_L3.fasta"
    if "advtig-reom1" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REO1_M1.fasta"
    if "advtig-reom2" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REO1_M2.fasta"
    if "advtig-reom3" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REO1_M3.fasta"
    if "advtig-reos1" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REO1_S1.fasta"
    if "advtig-reos2" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REO1_S2.fasta"
    if "advtig-reos3" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REO1_S3.fasta"
    if "advtig-reos4" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/REO1_S4.fasta"

    if "advtig-smrv" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/NC_001514.1.fna"
    if "advtig-betacoronvir" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/NC_006213.1.fna"
    if "advtig-mvm" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/NC_001510.1.fna"
    if "advtig-perv" == input_type.lower():
        spp_location = "/mnt/usersData/ADVTIG_fa/AF038600.1.fna"

    if "candida-18s-28s" == input_type.lower():
        spp_location = "/mnt/usersData/operons/library/EUK/NC_032096.1_Candida_albicans_SC5314_18S-ITS-28S_complement.fasta"
    if "aspergillus" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/any_fungus/library/fungi/GCF_000002855.3_ASM285v2_genomic_dustmasked.fna"
    if "malassezia" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/any_fungus/library/fungi/GCF_003290485.1_ASM329048v1_genomic_dustmasked.fna"
    if "malassezia globosa" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/any_fungus/library/fungi/GCF_000181695.1_ASM18169v1_genomic_dustmasked.fna"
    if "candida-albicans" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/any_fungus/library/fungi/GCF_000182965.3_ASM18296v3_genomic_dustmasked.fna"

    if "18s-aspergillus" == input_type.lower():
        spp_location = "/mnt/usersData/SILVA138/library/SILVA138/D63697.1.1733_Aspergillus_niger.fasta"
    if "18s-malassezia" == input_type.lower():
        spp_location = "/mnt/usersData/SILVA138/library/SILVA138/AF530542.1.1681_uncultured_Basidiomycota.fasta"
    if "18s-malassezia-restricta" == input_type.lower():
        spp_location = "/mnt/usersData/SILVA138/library/Clean_SSU_SILVA138/EU192367.1.1549_Malassezia_restricta.fasta"
    if "18s-candida-albicans" == input_type.lower():
        spp_location = "/mnt/usersData/SILVA138/library/Clean_SSU_SILVA138/AB013586.1.1769_Candida_albicans.fasta"

    if "cpn60-candida-albicans" == input_type.lower():
        spp_location = "/mnt/usersData/operons/library/cpn60/b7324_AY837556_Candida_albicans_ATCC_10231.fasta"

    if "16s-ec" == input_type.lower():
        spp_location = "/mnt/usersData/SILVA138/library/Clean_SSU_SILVA138/AB680517.1.1467_Escherichia_coli.fasta"
    if "16s-sa" == input_type.lower():
        spp_location = "/mnt/usersData/SILVA138/library/Clean_SSU_SILVA138/AB680000.1.1477_Staphylococcus_aureus.fasta"
    if "16s-pa7"== input_type.lower():
        spp_location = "/mnt/usersData/operons/library/BAC/NC_009656.1_Pseudomonas_aeruginosa_PA7,_16S-ITS-23S_complement.fasta"
    # if "16s-pa7"== input_type.lower():
    #     spp_location = "/mnt/usersData/SILVA138/library/Clean_SSU_SILVA138/"

    if "16s-cs" == input_type.lower():
        spp_location = "/mnt/usersData/SILVA138/library/Clean_SSU_SILVA138/KU950269.1.1399_Clostridium_sporogenes.fasta"
    if "clostridium-botulinum" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library/bacteria/GCF_000829015.1_ASM82901v1_genomic_dustmasked.fna"

    # FELV
    # spp_location = "/home/james/SequencingData/Centrifuge_libraries/viral/library/refseq_viral/GCF_000850105.1_ViralProj14686_genomic_dustmasked.fna"
    if "feline virus" == input_type.lower() or "felv" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/viral/FELV-RVDB.fasta"
    # PCV1
    if "porcine" == input_type.lower() or "pcv" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/viral/library/refseq_viral/GCF_000837765.1_ViralProj14053_genomic_dustmasked.fna" # PCV-1
    # E.coli K-12
    if "ecoli" == input_type.lower() or "escherichia" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library/bacteria/GCF_000974825.1_ASM97482v1_genomic_dustmasked.fna"
    # # MVM
    if "minute" == input_type.lower() or "mvm" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/viral/library/refseq_viral/GCF_000838465.1_ViralProj14019_genomic_dustmasked.fna"
    # felis catus
    if "cat" == input_type.lower() or "feline" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/cat/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna"
    # full human sequence
    if "human" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/human/human_input-sequences.fna"
    # full CHO sequence
    if "cho" == input_type.lower() or "chinese hamster" == input_type.lower():
        spp_location = "/mnt/usersData/krakenDB/library/vertebrate_mammalian/Chromosome/Cricetulus_griseus_strain_17A_GY-tax10029-GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"
    # full cutibacterium acnes sequence
    if "cacnes" == input_type.lower() or "cutibacterium acnes" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library/bacteria/GCF_000008345.1_ASM834v1_genomic_dustmasked.fna"
    # full staphylococcus epidermidis sequence
    if "sepidermidis" == input_type.lower() or "staphylococcus epidermidis" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library/bacteria/GCF_006742205.1_ASM674220v1_genomic_dustmasked.fna"
    if "staphylococcus aureus ca15" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library/bacteria/GCF_001021895.1_ASM102189v1_genomic_dustmasked.fna"
    if "saureus" == input_type.lower() or "staphylococcus aureus" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library/bacteria/GCF_003573855.1_ASM357385v1_genomic_dustmasked.fna"
    if "pao1" == input_type.lower() or "pseudomonas aeruginosa pao1" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library/bacteria/GCF_000006765.1_ASM676v1_genomic_dustmasked.fna"
    if "paeruginosa" == input_type.lower() or "pseudomonas aeruginosa" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library/bacteria/GCF_002075065.1_ASM207506v1_genomic_dustmasked.fna"
    if "lambda phage" == input_type.lower() or "lambda" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/ecoli-lambda-phage.fasta"
    if "klebsiella pneumoniae" == input_type.lower() or "kp" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library/bacteria/GCF_000009885.1_ASM988v1_genomic_dustmasked.fna"
    if "moraxella catarrhalis" == input_type.lower() or "mx" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library//bacteria_small_smart/GCF_000092265.1_ASM9226v1_genomic_dustmasked.fna"
    if "sordaria macrospora" == input_type.lower() or "sm" == input_type.lower():
        spp_location = "/home/james/SequencingData/NCBI_fungal_full/genomic/GCF_000182805.2_ASM18280v2_genomic.fna"
    if "lelliottia amnigena" == input_type.lower() or "la" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library//bacteria_small_smart/GCF_001652505.2_ASM165250v2_genomic_dustmasked.fna"
    if "bacillus subtilis" == input_type.lower() or "bs" == input_type.lower():
        spp_location = "/home/james/SequencingData/Centrifuge_libraries/bacteria/library//bacteria_small_smart/GCF_000146565.1_ASM14656v1_genomic_dustmasked.fna"

    bespoke_analyses = f"{directory}/analysis/bespoke/"
    os.makedirs(bespoke_analyses, exist_ok=True)

    samples = [line.rstrip("\n").split(",") for line in open(f"{sample_names}")]
    samples = [item for sublist in samples for item in sublist]

    clean_up = False
    no_split = True  # should be True when init-reads is True
    init_reads = False
    # processes = 5
    read_max_len = 150

    print(f"RNA status: {rna}")
    print(f"asm5 status: {asm5}")
    print(f"primary status: {primary}")

    with open(spp_location, "r") as f:
        first_line = f.readline()
        for_save = first_line.replace("\n", "")
        us_count = first_line.count('_')
        if us_count > 1:
            first_line = first_line.replace("_", " ")
        species_name = "-".join(first_line.split(" ")[1:])
        species_name = species_name.replace(",", "")
        species_name = species_name.replace("$", "")
        species_name = species_name.replace("\n", "")
        species_name = species_name.replace("/", "-")
        species_name = species_name.replace(".", "-")
        species_name = species_name.replace("(", "-")
        species_name = species_name.replace(")", "-")
        species_name = species_name.replace("[", "-")
        species_name = species_name.replace("]", "-")

    print(f"\nSpecies name: {species_name}\nLocation: {spp_location}\n")

    # no split will divide up the fasta into subsequences to zoom in on a location
    # use this for whole eukaroytic genomes like human or chinese hamster
    if no_split:
        main(
            samples,
            spp_location,
            directory,
            clean_up,
            species_name,
            bespoke_analyses,
            processes,
            read_max_len,
            init_reads,
            rna,
            asm5, primary,
            hum_fq,
            other
        )

    # this can be used when you want to look deeper at a specific region
    # chunking is relatively random however.
    if not no_split:
        temp_species_name = f"temp-{species_name}"
        temp_out = spp_location.split("/")[-1]
        new_name = for_save.split(" ")[0].split(".")[0]
        with tempfile.TemporaryDirectory() as path:
            # n = 1000
            n = 20
            no_lines = file_len(spp_location)
            print(no_lines)
            n_slices = round(no_lines / n)

            # check actual length of genome ... currently there is a discrepancy
            if os.stat(spp_location).st_size < 5e5:
                with open(spp_location, "r") as f:
                    fna = f.readlines()
                seq_only = fna[1:]
                seq_only = [line.replace("\n", "") for line in seq_only]
                length = 0
                for item in seq_only:
                    length += len(item)
                print(length)

            indices_to_add = list(range(n_slices, (n + 1) * n_slices, n_slices))
            additional_idx = round((no_lines - indices_to_add[-1]) / n_slices)

            for idx in range(additional_idx + 1):
                indices_to_add.append(indices_to_add[-1] + n_slices)

            new_spp_location = f"{path}/temp-{temp_out}"
            with open(spp_location, "r") as f_in, open(new_spp_location, "w") as f_out:
                for i, line in enumerate(f_in):
                    new_name_id = f"{new_name}{i}.1"
                    if i in indices_to_add:
                        f_out.write(f"{new_name_id} {species_name} {i}\n")
                        # print(f"{new_name_id} {species_name} {i}\n")
                    f_out.write(f"{line}")

            main(
                samples,
                new_spp_location,
                directory,
                clean_up,
                temp_species_name,
                bespoke_analyses,
                processes,
                read_max_len,
                init_reads,
                rna,
                asm5, primary,
                hum_fq,
                other
            )



# time ./pipeline/map_genome_to_reads.py -d "/mnt/usersData/Viral_human/" -s "/home/james/SMART-CAMP/configs/viral_DNA_all10.txt" -t PCV
# time ./pipeline/map_genome_to_reads.py -d "/mnt/usersData/Viral_human/" -s "/home/james/SMART-CAMP/configs/viral_DNA_all10.txt" -t FELV
# time ./pipeline/map_genome_to_reads.py -d "/mnt/usersData/Viral_human/" -s "/home/james/SMART-CAMP/configs/viral_DNA_all10.txt" -t human
# time ./pipeline/map_genome_to_reads.py -d "/mnt/usersData/Viral_CHO/" -s "/home/james/SMART-CAMP/configs/viral_DNA_all2.txt" -t MVM
# time ./pipeline/map_genome_to_reads.py -d "/mnt/usersData/Viral_CHO/" -s "/home/james/SMART-CAMP/configs/viral_DNA_all2.txt" -t CHO
