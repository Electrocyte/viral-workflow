#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 15:14:19 2021

@author: mangi
"""

encoding = 'utf-8'

from pathlib import Path
from shutil import copyfile
import json
import subprocess
import numpy as np
import argparse
import glob
import time
import os
import sys
import re
import pandas as pd
import shutil
import gzip
from typing import List
from functools import partial
from pipeline import Convert_U2T_fastq, mp_host_remove, \
     mp_trim_reads_v2, mp_demux_reads, database_config, \
     rank_centrifuge_predictions_for_aa_v2, \
     sample_name_sanity_check, \
     dep_flat_sample_classifier_test, mp_comb_ce_BN_v3, mp_cent_interpreter_v3, \
     dep_run_nanostat_analyses, guppy_demux
from q_control import adjust_q_scores, filter_q


def copy_and_rename(input_directory: str, s: str) -> None:
    # Search for files that match the specified pattern
    file_paths = glob.glob(f"{input_directory}/analysis/sample_data/{s}/trimmed/trimmed_{s}.fastq")
    if len(file_paths) > 0:
        file_path = file_paths[0]
        new_file_path = f"{input_directory}/analysis/sample_data/{s}/trimmed/original.fastq"

        # Copy the files to the new file names
        if not os.path.isfile(new_file_path):
            shutil.copy(file_path, new_file_path)


def rename_final_file(input_directory: str, s: str, quality: str, adjust_qscore: str) -> None:
    # Determine the name of the last file generated
    if adjust_qscore is not None and quality is not None:
        old_file_path = f"{input_directory}/analysis/sample_data/{s}/trimmed/t-adj-q{quality}.fastq"
    elif adjust_qscore is not None:
        old_file_path = f"{input_directory}/analysis/sample_data/{s}/trimmed/t-adjusted-{adjust_qscore}Q.fastq"
    elif quality is not None:
        old_file_path = f"{input_directory}/analysis/sample_data/{s}/trimmed/q{quality}.fastq"

    # New file path
    new_file_path = f"{input_directory}/analysis/sample_data/{s}/trimmed/trimmed_{s}.fastq"

    # copy the file
    shutil.copy(old_file_path, new_file_path)


def q_filters(samples: List, directory: str, adjust_qscore: str, q: str) -> None:
    # If both q is None and adjust_qscore is None, skip the operation
    if q is None and adjust_qscore is None:
        return

    # copy existing sample fastq files and rename (trimmed, demux, original)
    for s in samples:
        copy_and_rename(directory, s)

    # adjust down 1 q value each base
    if adjust_qscore is not None:
        adjust_q_scores.process_samples(directory, samples, adjust_qscore)

    if q is not None:
        filter_q.process_files(directory, q, adjust_qscore, samples)

    # renaming the final generated file to trimmed_{s}.fastq
    for s in samples:
        rename_final_file(directory, s, q, adjust_qscore)


# handle individual cat merged fastq files
def move_lonely_fastq(directory: str, sample_filename: str, fastq_file: str) -> str:
    if os.path.isfile(fastq_file):
        find_runid = subprocess.run(["rg", "--max-count", "1", '-i', "runid=", fastq_file], capture_output=True)
        stdout = str(find_runid.stdout, encoding)
        first_runid = re.compile(r"runid=(\w+)")
        is_runid = first_runid.findall(stdout)
        run_id = is_runid[0]

        os.makedirs(f'{directory}/{sample_filename}/', exist_ok=True)
        os.makedirs(f'{directory}/{sample_filename}/{run_id}/', exist_ok=True)
        os.makedirs(f'{directory}/{sample_filename}/{run_id}/fastq_pass/', exist_ok=True)

        destination = f'{directory}/{sample_filename}/{run_id}/fastq_pass/{sample_filename}.fastq'
        to_remove = f"{directory}/{sample_filename}"
        print(f"Moving file to: {destination}")
        os.rename(fastq_file, destination)
        return to_remove


# takes in pd.df with at least 'date', 'NA', 'strain', 'concentration_CFU', 'batch', 'duration_h'
def get_barcodes(df: pd.DataFrame) -> (list, list):
    print(df)
    df['Barcode'] = df['Barcode'].str.replace(" ", "")

    dates = df['date'].to_list()
    NAs = df['NA'].to_list()
    strains = df['strain'].to_list()
    CFUs = df['concentration_CFU'].to_list()
    batches = df['batch'].to_list()
    time = df['duration_h'].to_list()
    sample_id = [list(e) for e in zip(dates, NAs, strains, CFUs, batches, time)]
    identities = df['Identifier'].to_list()
    barcodes = df['Barcode'].to_list()

    final_samples = []
    for sample in sample_id:
        intermed_sample = []
        for item in sample:
            intermed_sample.append(str(item))
        intermed_sample = '_'.join(intermed_sample)
        final_samples.append(intermed_sample)

    folder_builder = [list(e) for e in zip(final_samples, identities, barcodes)]
    return final_samples, folder_builder

# split the barcoded samples by original identifier and demultiplex
# generate file directory for barcoded sample original names
# collect a list of intermediate folders to remove after analysis
def handle_barcoded_samples(directory: str, samples: list, builder: list, df: pd.DataFrame, threads: str, meta_fastq_fail: bool) -> (list, set):
    # find unique identifiers for all samples in list
    identifiers = []
    for index, sample in enumerate(samples):
        if builder[index][0] == sample:
            identifier = builder[index][1]
            identifiers.append(identifier)

    experiments = set(identifiers)

    available_experiments = []
    for check_experiment in experiments:
        if os.path.isdir (f'{directory}/{check_experiment}/'):
            available_experiments.append(check_experiment)
    available_experiments = set(available_experiments)

    print(samples, experiments, available_experiments)
    # build directory structure with identifiers
    list_tuples = []
    for experiment in available_experiments:
        if "Empty_DataFrame" not in str(directory):
            os.makedirs(f'{directory}/analysis/', exist_ok=True)
            os.makedirs(f'{directory}/analysis/sample_data/', exist_ok=True)
            os.makedirs(f'{directory}/analysis/sample_data/{experiment}/', exist_ok=True)
            os.makedirs(f'{directory}/analysis/sample_data/{experiment}/demux/', exist_ok=True)
            print(f'\n{directory}/{experiment}/', os.path.isdir (f'{directory}/{experiment}/'))

            # run qcat
            sub_experiment_df = df.loc[df.Identifier == experiment]
            kit = sub_experiment_df.Kit.unique()[0]

            barcode_file = f'{directory}/analysis/sample_data/{experiment}/demux/barcode*.fastq'
            print(f"Barcode: {barcode_file}")

            # check if demux files already exist or not
            # problem here is that when the file does not exist
            # glob returns an empty list, thus we check for this
            # and return an empty string if barcoded file
            # does not yet exist
            check_demux_complete = glob.glob(barcode_file)
            if len(check_demux_complete) > 0:
                check_demux_complete = check_demux_complete[0]
            if len(check_demux_complete) == 0:
                check_demux_complete = ""

            print(f"Barcode files: {check_demux_complete}")
            if not os.path.isfile (check_demux_complete):
                mp_demux_reads.main(threads, directory, experiment, meta_fastq_fail, kit)

            barcoded_files = glob.glob(f'{directory}/analysis/sample_data/{experiment}/demux/*.fastq')
            barcodes = [barcode.split('/')[-1].split('.')[0] for barcode in barcoded_files]
            for barcode in barcodes:
                if "none" in barcode:
                    barcodes.remove("none")
            list_tuples.append(tuple((barcodes, experiment)))

    # check the barcodes and identifiers against the original data table i.e. df
    removal = []
    for tuple_ in list_tuples:
        for barcode in tuple_[0]:
            sample_item = barcode_df[barcode_df['Barcode'].str.contains(barcode) & barcode_df['Identifier'].str.contains(tuple_[1])][['date', 'NA', 'strain', 'concentration_CFU', 'batch', 'duration_h']]
            dirty_str = sample_item.to_string(header=False, index=False, index_names=False).split('\n')
            clean_str = ['_'.join(ele.split()) for ele in dirty_str][0]
            if clean_str != "Empty_DataFrame":
                print(f'Generating files here: {directory}/analysis/sample_data/{clean_str}/')

                input_fastq_for_copying = f"{directory}/analysis/sample_data/{tuple_[1]}/demux/{barcode}.fastq"

                # use barcode copy for the new directory generation (will be removed after analysis)
                input_fastq_for_trimming = f"{directory}/analysis/sample_data/{tuple_[1]}/demux/{barcode}_copy.fastq"
                copyfile(input_fastq_for_copying, input_fastq_for_trimming)

                to_delete = move_lonely_fastq(directory, clean_str, input_fastq_for_trimming)
                removal.append(to_delete)

    print("\n\n\nFinishing barcoding")
    return removal, available_experiments


# generate file directory for analysis
def generate_file_dir(directory: str, sample: str) -> None:
    os.makedirs(f'{directory}/temp/', exist_ok=True)
    os.makedirs(f'{directory}/analysis/', exist_ok=True)
    os.makedirs(f'{directory}/analysis/sample_data/', exist_ok=True)
    os.makedirs(f'{directory}/analysis/sample_data/{sample}/', exist_ok=True)
    os.makedirs(f'{directory}/analysis/sample_data/{sample}/trimmed/', exist_ok=True)
    os.makedirs(f'{directory}/analysis/sample_data/{sample}/kraken/', exist_ok=True)
    os.makedirs(f'{directory}/analysis/sample_data/{sample}/centrifuge/', exist_ok=True)
    os.makedirs(f'{directory}/analysis/sample_data/{sample}/blastN/', exist_ok=True)


# demultiplex and trim reads
def run_read_trimming(samples: list, epoch_time: int, basecalled_already: bool, barcoded_sample: bool, barcoded_fastqs: str, top_dir: str, threads_: str, meta_fastq_fail: bool) -> None:
    no_samples = len(samples)
    n = 0
    for sample in samples:
        if not os.path.isdir(f'{top_dir}/analysis/sample_data/{sample}/blastN/'):
            generate_file_dir(top_dir, sample)
        if not os.path.isfile(f"{top_dir}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"):
            print(f"\n>>>Analysing {sample}<<<")
            n += 1
            print(f"Samples processed: {n}/{no_samples}")
            print(f"Sample: {sample}, is the sample basecalled: {basecalled_already}, barcoded_sample: {barcoded_sample}")

            if basecalled_already:
                print("\nMake a folder called 'temp' and add any basecalled fastq files here.\n")
                file_temp = f"{top_dir}/temp/{sample}.fastq"
                print(f"File to be moved {file_temp}")
                # for deletion is an optional choice used in the barcoding process, ignore it here in single sample processing
                move_lonely_fastq(top_dir, sample, file_temp)

            if Path(f'{top_dir}/{sample}').is_dir():
                # trim reads FAST.
                # step 1 - prepare chunked files if non-barcoded
                print("Running main one.")
                mp_trim_reads_v2.main(threads_, top_dir, sample, meta_fastq_fail, barcoded_sample, epoch_time)

    # step 2 - run porechop on chunked files
    if not barcoded_sample:
        print("Running main two.")
        mp_trim_reads_v2.main_two(threads_, samples, top_dir, meta_fastq_fail)


def run_host_removal(q: str, adjust_qscore: str, top_dir: str, samples: list, hosts: list, lenhosts: int, processes_input: str, cent_index: str, hsbn_path_str: str, blastN_lib, krakenuniq_command: str, github: str, epoch_dir_host: str, kingdom: str, blastN_libraries: dict) -> (bool, list, bool, str):
# complete conversion of U to T, remove host reads and commence metagenomic classification
    new_fqs = []
    host_bool = False
    no_host = ""
    Direct_RNA = False
    new_fq = ""
    ce_host_read_count = 0
    hsbn_host_read_count = 0
    kr_host_read_count = 0

    for sample in samples:
        print(f"\n\n\nCurrent sample: {sample}")
        if Path(f'{top_dir}/analysis/sample_data/{sample}').is_dir():
            if Path(f'{top_dir}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq').is_file():

                trim_dir = f'{top_dir}/analysis/sample_data/{sample}/trimmed/'
                trim_fq = f'{trim_dir}/trimmed_{sample}.fastq'

                if "-RNA-" in sample or "_RNA_" in sample:
                    Convert_U2T_fastq.main(trim_fq, trim_dir, sample)

                if lenhosts > 0:
                    no_host = "_no_host"
                    for host in hosts:
#                        if host in sample:
                        host_bool = True
                        if host in "CHO":
                            host = "CHO"
                        elif host in "barramundi":
                            host = "barramundi"
                        else:
                            host = "human"
                        print(f"Host found: {host}")
                        print(f"Removing host reads from sample. Example: {top_dir}/analysis/sample_data/{sample}/trimmed/no_host_{sample}.fastq")
                        new_fq, ce_host_read_count, hsbn_host_read_count, kr_host_read_count = mp_host_remove.remove_host_reads(sample, top_dir, processes_input, cent_index, Direct_RNA, hsbn_path_str, host, epoch_dir_host, github, krakenuniq_command, blastN_libraries)
                        new_fqs.append(new_fq)
                        new_fq = ""

                kr_reads_list = [int(kr_host_read_count)]
                with open(f"{epoch_dir_host}/kraken-host-reads-{sample}.json", 'w', encoding='utf-8') as f:
                    json.dump(kr_reads_list, f)
                BLAST_reads_list = [int(hsbn_host_read_count)]
                with open(f"{epoch_dir_host}/BLAST-host-reads-{sample}.json", 'w', encoding='utf-8') as f:
                    json.dump(BLAST_reads_list, f)
                Ce_reads_list = [int(ce_host_read_count)]
                with open(f"{epoch_dir_host}/Centrifuge-host-reads-{sample}.json", 'w', encoding='utf-8') as f:
                    json.dump(Ce_reads_list, f)
    print(f"Files for metagenomic classification: {new_fqs}")
    return host_bool, new_fqs, Direct_RNA, no_host


def run_metagenome_classifiers(top_dir: str, hsbn_path_str: str, clf_idx: str, bN_name: str, blastN_lib: str, \
                               processes_input: int, kingdom: str, no_host, new_fqs: str, \
                               epoch_dir_host: str, github: str, samples: list, skip: bool) -> (list):

    print(f"\n\n\n\n{clf_idx}{kingdom}\n\n")
    for i in range(2):
        print(f"{clf_idx}{kingdom}")
    dep_flat_sample_classifier_test.main(top_dir, hsbn_path_str, clf_idx, bN_name, blastN_lib, processes_input, kingdom, no_host, new_fqs, epoch_dir_host, github, samples, skip)

    host_background_BLAST = 0
    ce_classified_reads = 0
    hs_classified_reads = 0
    hs_reads = 0
    ce_reads = 0
    hsbn_host_read_count = 0
    ce_host_read_count = 0
    total_reads = []
    for sample in samples:

        # import BLAST data
        if not skip:
            with open(f"{epoch_dir_host}/BLAST-AA-{sample}.json") as json_b_data:
                BLAST_json = json.loads(json_b_data.read())
            hs_reads, hs_classified_reads = BLAST_json # (int, int)

            # check if host reads present:
            if len(no_host) > 0:
                host_BLAST_file = f"{epoch_dir_host}/BLAST-host-reads-{sample}.json"
                if os.path.isfile(host_BLAST_file):
                    with open(host_BLAST_file) as json_b_host_data:
                        BLAST_host_json = json.loads(json_b_host_data.read())
                    hsbn_host_read_count = BLAST_host_json[0] # (int)
                    hsbn_host_read_count = int(hsbn_host_read_count)
        # import centrifuge data
        with open(f"{epoch_dir_host}/Centrifuge-AA-{sample}.json") as json_c_data:
            Centrifuge_json = json.loads(json_c_data.read())
        ce_reads, ce_classified_reads = Centrifuge_json # (int, int)

        # check if host reads present:
        if len(no_host) > 0:
            host_centrifuge_file = f"{epoch_dir_host}/Centrifuge-host-reads-{sample}.json"
            if os.path.isfile(host_centrifuge_file):
                with open(host_centrifuge_file) as json_c_host_data:
                    Centrifuge_host_json = json.loads(json_c_host_data.read())
                ce_host_read_count = Centrifuge_host_json[0] # (int)
                ce_host_read_count = int(ce_host_read_count)

        if hs_reads > 0:
            host_background_BLAST = (hsbn_host_read_count / hs_reads * 100)
        total_reads.append([sample, hs_reads, hs_classified_reads, ce_reads, ce_classified_reads, hsbn_host_read_count, ce_host_read_count, host_background_BLAST])
        print(f"\nHost reads: Centrifuge: {ce_host_read_count}, BLAST: {hsbn_host_read_count}, all BLAST: {hs_reads}, classified BLAST: {hs_classified_reads}, all Centrifuge: {ce_reads}, classified Centrifuge: {ce_classified_reads}")
        if not Path(f'{top_dir}/{sample}').is_dir():
            samples.remove(sample)
    return total_reads



def combine_cent_BLAST_wrapper(directory: str, meta_data_name: str, barcoded: bool, kingdom: str, \
                                bN_name: str, no_host: str, targ_species: dict, \
                                blastN_lib: dict, samples: list, skip: bool) -> (str, str, list, list):

    cent_df_str = mp_cent_interpreter_v3.main(directory, kingdom, no_host, targ_species)

    print(f"Centrifuge save file: {cent_df_str}")
    # analyse the ranks for the centrifuge analysis.

    # handle a single adventitious agent test species
    cent_spp_ranked = []
    analysis_directory = f"{directory}/analysis/"

    if not skip:

        updated_meta_data_name, agg_df_name, BLAST_spp_ranked = mp_comb_ce_BN_v3.main(directory, barcoded, kingdom, bN_name, meta_data_name, no_host, cent_df_str, targ_species)

        print(f"New meta file: {updated_meta_data_name}, analysis type: {agg_df_name}")

        # return updated_meta_data_name, agg_df_name, BLAST_spp_ranked, cent_spp_ranked


def de_gunzip_file(files: list) -> None:
    for file in files:
        with gzip.open(file, 'rb') as f_in:
            not_gz = file.replace(".gz", "")
            if not os.path.isfile(not_gz):
                with open(not_gz, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        print(file, file.replace(".gz", ""))


def move_files_based_on_criteria(src_directory, dst_directory):
    # Ensure only FASTQ files (.fastq or .fastq.gz) are moved
    for filename in os.listdir(src_directory):
        if filename.endswith('.fastq') or filename.endswith('.fastq.gz'):
            src_file = os.path.join(src_directory, filename)
            dst_file = os.path.join(dst_directory, filename)
            if os.path.isfile(src_file):  # Ensure it's a file, not a directory
                shutil.move(src_file, dst_file)


parser = argparse.ArgumentParser(description='Running FAST metagenomic assessor')

# input file dir
parser.add_argument('-d', '--directory', help='Input file directory to reference genome checked FASTQ files.')

parser.add_argument('-q', '--qscore', default=None, help='Filter on qscore.')

# input file dir
parser.add_argument('-ts', '--target_species', help='''Name of species per hs-blastn index. This will be used to build a machine learning model that will reduce the number of potential predictions.
                    Only takes, dictionary, including a single species e.g. '{"Pseudomonas aeruginosa":"PA","Cutibacterium acnes":"Pacnes","Escherichia coli":"EC"}'.
                    The key is the species and genus, the value is the sample label including the substring that should be common to all your samples.
                    In the example given, samples were spiked with Pseudomonas aeruginosa and the sample name was PAO1unfiltered, thus PA will be a substring.''')

# input centrifuge index location
# might be discontinued
parser.add_argument('-ci', '--centrifuge_index', help='Centrifuge classifier index directory location e.g. ~/SequencingData/Centrifuge_libraries/')

# Number of threads to use \
parser.add_argument('-t', '--threads', help='Number of threads / CPU cores to deploy. Default 4.', default=4, type=int)

# # Quality score
# parser.add_argument('-q', '--q_score', help='Read quality score. Default 7.', default=7, type=int)

# plot time to next read
parser.add_argument('-r', '--next_read', action='store_true', help='Time to next read.')

# plot blastN output
parser.add_argument('-p', '--plots', action='store_true', help='Plots for blastN analysis (or not).')

# host species
parser.add_argument('-hn', '--host_name', help="Hosts in sample e.g. TC for T cell or WBC for white blood cell, in format 'WBC,TC,Jurkat,CHO' ", default='')

# include fastq fail files
parser.add_argument('-qf', '--fastq_fail', action='store_true', help="Include fastq fail files defined by having quality score <7; -qf if true")

# include fastq fail files
parser.add_argument('-a', '--adjust_qscore', default=None, help="Adjust q-score by n")

# meta data file
parser.add_argument('-m', '--meta_data', action='store_true', help="Use priority candidates or not, set to -m if true")

# meta data file
parser.add_argument('-bin', '--bin_reads', action='store_true', help="Bin true and false positive reads, set to -m if true")

# meta data file
parser.add_argument('-skip', '--skip_main', action='store_true', help="Only do centrifuge metagenomic classification")

# include fastq fail files
parser.add_argument('-ha', '--host_analyses', action='store_true', help="Analyse host data using krakenuniq -- generates kmer coverage and percent coverage.")

# input file list of str
parser.add_argument('-fd', '--fastq_dirs', help='Input text file with sample names of interest to reference genome checked FASTQ files.')

# input file are barcoded
parser.add_argument('-bc', '--barcodes', help='Input sample file table with columns date [20210220], \nNA [nucleic acid type e.g. RNA], strain or sample [str], concentration_CFU [1000CFU], \nbatch [1], duration_h [int], Flowcell_no [if you want to know this], \nSample_original [not important], Identifier [str e.g. S5B1], \nBarcode [e.g. barcode01], Total_cycles [PCR cycles e.g. 0, 36], Spike_species [e.g. Pseudomonas]')

random_spp_ = \
"""If you wish to improve processing speed, you can randomly sample from the analysed refseq species.
Default is set to 1 refseq files for a given identified species.
This can drastically enhance performance when some species, such as E.coli, can have 1000s of refseq genomes."""

# if there are many refseq samples, default is to choose 10 randomly
parser.add_argument('-rs', '--random_spp', help=random_spp_, default=1, type=int)

# clean up after running entire script by deleting intermediate files
parser.add_argument('-c', '--clean', action='store_true', help='Clean up after running entire script by deleting intermediate files. Set to -c if True.')

# pre-basecalled reads
parser.add_argument('-b', '--basecalled', action='store_true', help='Reads already basecalled? Set to -b if True.')

############### GLOBAL VARIABLES ###############

blastN_libraries = database_config.blastN_databases()
coverage_libraries = database_config.coverage_databases()
host_coverage_libraries = database_config.host_coverage_databases()

# choose a blastN library
parser.add_argument('-bl', '--blastN_lib', choices=blastN_libraries.keys(), default="16S_16S")#, help='Input desired library for blastN analysis e.g. 16S_16S, 16S_23S, NR99_SSU, NR99_LSU, SILVA, bacteria, human, cviral, virus, uviral, fungal_all, small_bacteria, filter_bacteria, chinese_hamster', default="16S_16S")

######### !!!!!!!!!!!!!!!!!!!! #########
# this may require changing if not calling the function from GitHub directory
# cwd = os.getcwd()
github = "/home/james/SMART-CAMP/"
krakenuniq_command = "/home/james/krakenuniq/krakenuniq/krakenuniq"
hsbn_path_str = "/home/james/apps/hs-blastn"
######### !!!!!!!!!!!!!!!!!!!! #########

# must parse args FIRST
args = parser.parse_args()

q = args.qscore

top_dir = args.directory
targ_species = args.target_species

dir_input = Path(f'{args.directory}')
cent_index = args.centrifuge_index
threads_ = str(args.threads)
# q_score_ = str(args.q_score)

# boolean inputs
########################
to_plot  = args.plots
host_analysis = args.host_analyses
bin_AA_reads = args.bin_reads
skip = args.skip_main
timetonextread = args.next_read
clean_up = args.clean
basecalled_already  = args.basecalled
meta_data_bool = args.meta_data
adjust_qscore = args.adjust_qscore
hosts = args.host_name
meta_fastq_fail = args.fastq_fail
########################

lenhosts = len(hosts)
hosts = hosts.split(",")

# what is the purpose of this check?
if meta_data_bool:
    meta_data_name = "priority"

bN_name = args.blastN_lib

blastN_lib = blastN_libraries[args.blastN_lib]
random_spp_sample = args.random_spp
print(f"Library for blastN metagenomic analysis: {blastN_lib}, bN_name: {bN_name}")

# filename with directory
multiple_fastq = args.fastq_dirs

# meta table containing barcode information, requires barcode, identifier, date, nucleic acid (NA), strain, CFU
barcoded_fastqs = args.barcodes

if not isinstance(top_dir, list):
    top_dir = Path(f'{top_dir}' )

############### GLOBAL VARIABLES ###############

folder_builder = []
must_be_removed = []
experiment_names = set()
barcoded_sample = False

if barcoded_fastqs != None:
    barcoded_sample = True
    print("Barcoded preprocessing")
    barcode_df = pd.read_csv(f'{barcoded_fastqs}')

    samples, folder_builder = get_barcodes(barcode_df)
    print(folder_builder)

    unique_barcode = "barcode144"
    if unique_barcode in barcode_df["Barcode"].values:
        print(unique_barcode in barcode_df["Barcode"].values)
        for fb in folder_builder:
            samp, expt, barc = fb
            pass_folders = glob.glob(f"{top_dir}/{expt}/**/fastq_pass/", recursive=True)
            fail_folders = glob.glob(f"{top_dir}/{expt}/**/fastq_fail/", recursive=True)
            print(f"\nPass folder: {pass_folders}, fail folder: {fail_folders}")

            for pass_folder in pass_folders:
                pass_unique_folder = os.path.join(pass_folder, unique_barcode)
                os.makedirs(pass_unique_folder, exist_ok=True)
                move_files_based_on_criteria(pass_folder, pass_unique_folder)

            for fail_folder in fail_folders:
                fail_unique_folder = os.path.join(fail_folder, unique_barcode)
                os.makedirs(fail_unique_folder, exist_ok=True)
                move_files_based_on_criteria(fail_folder, fail_unique_folder)

            print(f"Pass barcode folder: {pass_unique_folder}, fail barcode folder: {fail_unique_folder}\n")

    guppy_demux_idx = barcode_df.loc[barcode_df["Sample_original"].str.contains("guppy-demux")].index
    qcat_demux_samples_idx = list(set(list(barcode_df.index)) - set(guppy_demux_idx))
    qcat_demux_samples = barcode_df.iloc[qcat_demux_samples_idx]
    qcat_demux_samples["sample"] = qcat_demux_samples['date'].astype(str) \
    +"_"+ qcat_demux_samples['NA'] +"_"+ qcat_demux_samples['strain'] \
        +"_"+ qcat_demux_samples['concentration_CFU'].astype(str) \
            +"_"+ qcat_demux_samples['batch'].astype(str) \
                +"_"+ qcat_demux_samples['duration_h'].astype(str)
    qcat_demux_samples_check = list(qcat_demux_samples["sample"].unique())

    guppy_demux_samples = barcode_df.iloc[guppy_demux_idx]
    guppy_demux_samples["sample"] = guppy_demux_samples['date'].astype(str) \
        +"_"+ guppy_demux_samples['NA'] +"_"+ guppy_demux_samples['strain'] \
            +"_"+ guppy_demux_samples['concentration_CFU'].astype(str) \
                +"_"+ guppy_demux_samples['batch'].astype(str) \
                    +"_"+ guppy_demux_samples['duration_h'].astype(str)
    guppy_demux_sample_check = list(guppy_demux_samples["sample"].unique())
    folder_builder_guppy_demux = [i for n,i in enumerate(folder_builder) if folder_builder[n][0] in guppy_demux_sample_check]

    print("Running qcat demux")
    must_be_removed, experiment_names = handle_barcoded_samples(top_dir, qcat_demux_samples_check, folder_builder, barcode_df, threads_, meta_fastq_fail)

    if unique_barcode in barcode_df["Barcode"].values:
        print(unique_barcode in barcode_df["Barcode"].values)
        for fb in folder_builder:
            samp, expt, barc = fb
            out_ss_loc = f"{top_dir}/{samp}/{barc}/"
            ss_loc = f"{top_dir}/{expt}/**/seq*summ*t"
            print(ss_loc)
            os.makedirs(out_ss_loc, exist_ok=True)
            ss_file = glob.glob(ss_loc, recursive=True)
            ss_out = glob.glob(out_ss_loc, recursive=True)
            if ss_out:
                for so in ss_out:
                    print(so, ss_file)
                    if len(ss_file) > 0:
                        print(ss_file)
                        if os.path.isfile(ss_file[0]):  # Ensure it's a file, not a directory
                            print(f"Summary file: {ss_file[0]}, copying to: {so}")
                            shutil.copy(ss_file[0], so)
                            saved_summ = glob.glob(f'{so}/seq*summ*t')[0]
                            print(f"{os.stat(saved_summ).st_size/1e6:.0f} MB")

    dep_run_nanostat_analyses.run_nanostat(top_dir, threads_, barcoded_fastqs)

    prefix = '/'.join(barcoded_fastqs.split("/")[:-1])
    suffix = barcoded_fastqs.split("/")[-1].split(".")[0]
    temp_file = f"{prefix}/temp_{suffix}.txt"

    # use this for further processing regardless

    print(f"Generating temp file: {temp_file}")
    with open(temp_file, 'w') as f:
        for sample in samples:
            f.write(f"{sample}\n")

    print("Running guppy demux")
    guppy_demux.main(folder_builder_guppy_demux, top_dir)
    multiple_fastq = temp_file
    print(f"\n\nTemporary file name: {multiple_fastq}\n\n")
    if meta_fastq_fail:
        print(f'Please run: time /path/to/mp_metagenomic_assessment_v4.py -d {top_dir} -t {threads_} -ci {cent_index} -c -qf -bl {blastN_lib} -m -rs {random_spp_sample} -fd {multiple_fastq} -hn {hosts} ')
    else:
        print(f'Please run: time /path/to/mp_metagenomic_assessment_v4.py -d {top_dir} -t {threads_} -ci {cent_index} -c -bl {blastN_lib} -m -rs {random_spp_sample} -fd {multiple_fastq} -hn {hosts[0]} ')
    print("Exiting barcode run.")
    sys.exit()
        # dep_run_nanostat_analyses.run_nanostat(top_dir, threads_, temp_file)

if multiple_fastq != None:
    print("Multiple non-barcoded samples for preprocessing")
    samples = [line.rstrip('\n').split(',') for line in open(f'{multiple_fastq}')]
    samples = [item for sublist in samples for item in sublist]

    for sample in samples:
        gz_files = glob.glob(f"{top_dir}/{sample}/**/fastq_*/*q*z")
        if len(gz_files) > 0:
            de_gunzip_file(gz_files)

    dep_run_nanostat_analyses.run_nanostat(top_dir, threads_, multiple_fastq)

# check sample names adhere to the convention and don't break the rest of the analysis.
print(samples)
sample_name_sanity_check.check_sample_name(samples)

epoch_time = int(time.time())
epoch_dir = f"{top_dir}/analysis/run/{str(epoch_time)}"
epoch_dir_host = f"{epoch_dir}/host/"
os.makedirs(f"{top_dir}/analysis/run/", exist_ok=True)
os.makedirs(epoch_dir, exist_ok=True)
os.makedirs(epoch_dir_host, exist_ok=True)

print(f"\n{bN_name} {bN_name == 'uviral25'} {bN_name == 'uviral25-2'}\n")

kingdom = ""
if bN_name == "bacteria" or bN_name == "16S_16S" or bN_name == "16S_23S" or bN_name == "filter_bacteria":
    kingdom = "bacteria"
if bN_name == "c-viral29":
    kingdom = "c-viral29"
if bN_name == "uviral29":
    kingdom = "uv29"
if bN_name == "uviral25-2":
    kingdom = "mini-u-viral2"
if bN_name == "uviral25":
    kingdom = "mini-u-viral"
if bN_name == "virus":
    kingdom = "v_f_b2"
if bN_name == "cviral25":
    kingdom = "mini-c-viral2"
if bN_name == "fungal_all":
    kingdom = "fungus"
if bN_name == "v_f_b":
    kingdom = "v_f_b"
if bN_name == "v_f_b2":
    kingdom = "v_f_b2"
if bN_name == "v_f_b3":
    kingdom = "v_f_b3"
if bN_name == "v_f_b4":
    kingdom = "v_f_b4"
# if bN_name == "v_f_b2_watson_BiSu":
#     kingdom = "v_f_b2_watson_BiSu"

print(f"Library for centrifuge metagenomic analysis: {cent_index}, kingdom: {kingdom}")

latest_run = ""

# if not skip:
print("Trimming reads")
run_read_trimming(samples, epoch_time, basecalled_already, barcoded_sample, barcoded_fastqs, top_dir, threads_, meta_fastq_fail)

q_filters(samples, top_dir, adjust_qscore, q)

print("\nRunning metagenomic classification.")
host_bool, new_fqs, Direct_RNA, no_host = run_host_removal(q, adjust_qscore, top_dir, samples, hosts, lenhosts, threads_, cent_index, hsbn_path_str, blastN_lib, krakenuniq_command, github, epoch_dir_host, kingdom, blastN_libraries)

total_reads = run_metagenome_classifiers(top_dir, hsbn_path_str, cent_index, bN_name, blastN_lib, threads_, kingdom, no_host, new_fqs, epoch_dir_host, github, samples, skip)

if not bN_name == "uviral25" or not bN_name == "cviral25" or not bN_name == "uviral25-2" or not bN_name == "uviral29":

    combine_cent_BLAST_wrapper(top_dir, meta_data_name, barcoded_sample, kingdom, bN_name, no_host, targ_species, blastN_lib, samples, skip)

    nanodf = pd.read_csv(f"{top_dir}/analysis/nanoplot_summary_data.csv")
    print(nanodf)
