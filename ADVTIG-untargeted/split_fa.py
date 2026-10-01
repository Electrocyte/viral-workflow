#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from itertools import zip_longest
from typing import Iterable, Optional, Tuple
import pandas as pd
import sys
import re
import json


def grouper(iterable: Iterable, n: int, fillvalue: Optional[str] = None) -> Iterable[Tuple]:
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


# Paths to the input files
fasta_file = "/home/james/SequencingData/NCBI_RVDB/virus_RVDBs/C-RVDBv29.0-fixed4.fasta"
mapping_file = "/home/james/SequencingData/NCBI_RVDB/virus_RVDBs/fixed_CRVDB_V29_seqID.map"
output_dir = "/home/james/SequencingData/NCBI_RVDB/virus_RVDBs/CRVDB_V29_split_fa/"

fasta_file = "/mnt/usersData/NCBI_RVDB/virus_RVDBs/U-RVDBv29.0-fixed2.fasta"
mapping_file = "/mnt/usersData/NCBI_RVDB/virus_RVDBs/fixed_URVDB_V29_seqID.map"
output_dir = "/mnt/usersData/NCBI_RVDB/virus_RVDBs/URVDB_V29_split_fa/"

fasta_file = "/mnt/usersData/plsdb/plasmid-sequences-fixed.fna"
mapping_file = "/mnt/usersData/plsdb/plasmid_seqid2taxid.map"
output_dir = "/mnt/usersData/plsdb/plasmid_split_fa/"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load the mapping file into a DataFrame
print(f"Reading mapping file: {mapping_file}")
mapping_df = pd.read_csv(mapping_file)
print("Columns detected:", mapping_df.columns.tolist())

code_to_name = dict(zip(mapping_df["code"], mapping_df["name"]))

index_dict = {}

# Process the FASTA file and split into individual files
if os.path.isfile(fasta_file):
    with open(fasta_file, "r") as infile:
        for idx, (header, sequence) in enumerate(grouper(infile, 2, fillvalue=None)):
            if header is None or sequence is None:
                continue  # Skip incomplete pairs

            header = header.strip()  # Remove any trailing newline characters
            sequence = sequence.strip()

            if not header.startswith(">"):
                print(f"rg -i -A 2 -B 2 {header} {fasta_file}")
                raise ValueError(f"Invalid FASTA header: {header}")

            # original
            # seq_id = header[1:]  # Remove the '>' to extract the sequence ID
            # hopefully more agnostic
            seq_id = header[1:].split()[0]

            name = code_to_name.get(seq_id, seq_id)  # Get the descriptive name or fallback to seq_id
            # print(name)
            # print(seq_id)

            if seq_id == "acc":
                print(header)
                print(sequence)
            # Generate the output file path
            try:
                out_name = name.replace(' ', '_').replace('/', '_')[:40].replace(".","_")
                out_name = re.sub(r'[^\w\s_]', '', out_name)
                # print(out_name)
            except:
                print(header)
                print(sequence[:20])
                print(f"Error with {name}")
                sys.exit(1)
            save_name = f"{seq_id}_{out_name}.fasta"
            output_file = os.path.join(output_dir, save_name)

            index_dict[seq_id] = save_name

            # # Write the sequence to the output file
            with open(output_file, "w") as outfile:
                outfile.write(f"{header}\n{sequence}\n")

            print(f"Created file {output_file}")
            # break

print("All sequences have been split and saved.")

# save index to json
# json_save = "/home/james/SequencingData/NCBI_RVDB/virus_RVDBs/U-RVDB-clean-fix2-index.json"
json_save = "/home/james/SequencingData/NCBI_RVDB/virus_RVDBs/C-RVDB-V29-index.json"
with open(json_save, "w") as f:
    json.dump(index_dict, f)
print(f"Index saved to json: {json_save}")
