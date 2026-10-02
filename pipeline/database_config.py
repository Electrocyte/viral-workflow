#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 11:42:01 2021

@author: mangi
"""

from pathlib import Path
import os


# function that outputs dictionary of coverage databases, this can be added to if you custom databases are required.
def coverage_databases() -> dict:
    # human sequence
    chinese_hamster_seq = "/mnt/usersData/krakenDB/vertebrate/vertebrate_mammalian/Chromosome/chinese_hamster/"

    barramundi = "/mnt/usersData/NCBI_Barramundi/GCF_001640805.1_ASM164080v1_genomic.fna"

    # full human sequence
    full_human_seq = f"{str(Path.home())}/SequencingData/Centrifuge_libraries/human/"

    coverage_libraries = {
        "human": full_human_seq,
        "barramundi": barramundi,
        "chinese_hamster": chinese_hamster_seq,
    }
    return coverage_libraries


def blastN_databases() -> dict:
    # viral-bacterial-fungal SEQUENCES (uses NCBI virus, fungal and smart small bacteria)
    v_f_b = (
        f"{str(Path.home())}/SequencingData/Centrifuge_libraries/v_f_b/v-f-b-small.fna"
    )

    # viral-bacterial-fungal SEQUENCES (uses NCBI virus, fungal and smart small bacteria)
    # includes clostridium sporogens and candida albicans
    v_f_b2 = (
        f"{str(Path.home())}/SequencingData/Centrifuge_libraries/v_f_b/v-f-b-small2.fna"
    )

    v_f_b3 = (
        f"{str(Path.home())}/SequencingData/Centrifuge_libraries/v_f_b/v-f-b-small3.fna"
    )

    v_f_b4 = (
        f"{str(Path.home())}/SequencingData/Centrifuge_libraries/v_f_b/v-f-b-small4.fna"
    )

    v_f_b2_w_BiSu = (
        f"{str(Path.home())}/SequencingData/Centrifuge_libraries/v_f_b/v-f-b-small2-watson-BiSu.fna"
    )

    uv29 = (f"/mnt/usersData/NCBI_RVDB/virus_RVDBs/U-RVDBv29.0-fixed2.fasta")

    un_viral_seq = (
        f"/mnt/usersData/NCBI_RVDB/virus_RVDBs/U-RVDB-clean-fix.fasta"
    )

    un_viral_seq2 = (
        f"/mnt/usersData/NCBI_RVDB/virus_RVDBs/U-RVDB-clean-fix2.fasta"
    )

    c_viral_seq3 = (
        f"/mnt/usersData/NCBI_RVDB/virus_RVDBs/C-RVDBv29.0-fixed4.fasta"
    )

    # # clustered viral
    cl_viral_seq25 = (
        f"/home/james/SequencingData/Centrifuge_libraries/viral/cRVDB2-input-sequences.fna"
    )

    # CH sequence
    chinese_hamster_seq = "/mnt/usersData/krakenDB/vertebrate/vertebrate_mammalian/Chromosome/chinese_hamster/GCF_000223135.1_CHOK1_CriGri_1.0_genomic.fna"

    barramundi = "/mnt/usersData/NCBI_Barramundi/GCF_001640805.1_ASM164080v1_genomic.fna"

    # human sequence
    full_human_seq = (
        f"/mnt/usersData/NCBI_human/fasta/GRCh38_latest_genomic.fna"
    )


    print("Checking if high speed blastN files exist:")
    print(f"v_f_b: {os.path.isfile(v_f_b)}")
    print(f"v_f_b2: {os.path.isfile(v_f_b2)}")
    print(f"v_f_b3: {os.path.isfile(v_f_b3)}")
    print(f"v_f_b4: {os.path.isfile(v_f_b4)}")
    print(f"v_f_b2: {os.path.isfile(v_f_b2_w_BiSu)}")
    print(f"Human all: {os.path.isfile(full_human_seq)}")
    print(f"CHO: {os.path.isfile(chinese_hamster_seq)}\n")


    blastN_libraries = {
        "v_f_b": v_f_b,
        "v_f_b2": v_f_b2,
        "v_f_b3": v_f_b3,
        "v_f_b4": v_f_b4,
        "v_f_b2_watson_BiSu": v_f_b2_w_BiSu,
        "human": full_human_seq,
        "chinese_hamster": chinese_hamster_seq,
        "barramundi": barramundi,
        "uv29":uv29,
        "uviral25":un_viral_seq,
        "uviral25-2":un_viral_seq2,
        "c-viral29":c_viral_seq3,
        "cviral25": cl_viral_seq25,
    }
    return blastN_libraries


# function that outputs dictionary of coverage databases, this can be added to if you custom databases are required.
def host_coverage_databases() -> dict:

    # CHO clean up
    # find ONLY the chromosome headers
    # rg -i ">.+chromosome" Cricetulus_griseus_strain_17A_GY-tax10029-GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna
    # find all unplaced contigs
    # rg -i ">.+unplaced" Cricetulus_griseus_strain_17A_GY-tax10029-GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna

    # delete all lines after given line - these are unassigned scaffolds
    # sed '28715648,$d' Cricetulus_griseus_strain_17A_GY-tax10029-GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna > CHO_chromosomes_only.fna

    # delete lines after line number -- chromosome 1
    # sed '3446229,$d' CHO_chromosomes_only.fna > CHO_chromosomes_only_1A.fna
    # delete lines before line number -- all other chromosomes
    # sed '6876127,$!d' CHO_chromosomes_only.fna > CHO_chromosomes_only_1B.fna

    chinese_hamster_chromosomes = "/mnt/usersData/krakenDB/library/vertebrate_mammalian/Chromosome/Cricetulus_griseus_strain_17A_GY-tax10029-GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"

    barramundi = "/mnt/usersData/NCBI_Barramundi/GCF_001640805.1_ASM164080v1_genomic.fna"

    # full human sequence
    full_human_seq = "/home/james/SequencingData/Centrifuge_libraries/human/human_input-sequences.fna"

    print("Checking if coverage folders exist:")
    print(f"Full human: {Path(full_human_seq).is_file()}")
    print(
        f"Chinese_hamster_seq (chromosomes): {Path(chinese_hamster_chromosomes).is_file()}"
    )

    host_coverage_libraries = {
        "human": full_human_seq,
        "chinese_hamster_chromosomes": chinese_hamster_chromosomes,
        "barramundi": barramundi,
    }
    return host_coverage_libraries
