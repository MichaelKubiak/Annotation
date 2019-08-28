#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to identify Swissprot entries with EC numbers
# ------------------------------------------------------------------------------------------------------

from pathlib import Path
import re
import pandas as pd

# ------------------------------------------------------------------------------------------------------
# read in files

with open(str(Path.home()) + "/Annotation/Data/enzyme.dat") as efile:
    enzyme = efile.readlines()

with open(str(Path.home()) + "/Annotation/Data/uniprot_sprot.fasta") as sprotfile:
    sprot = sprotfile.readlines()

# ------------------------------------------------------------------------------------------------------
# get a dictionary of EC numbers and the IDs of proteins of that classification

ECs = {}
current = ""
for line in enzyme:
    if line.startswith("ID"):
        current = line.split()[1]
        ECs[current] = []
    elif line.startswith("DR"):
        for protein in re.finditer(r"\w+,", line):

            ECs[current].append(protein.group().split(",")[0])

# ------------------------------------------------------------------------------------------------------
# get a dictionary of protein IDs with their sequences

proteins = {}
current = ""
for line in sprot:
    line = line.rstrip("\n")
    if line.startswith(">"):
        current = line.split("|")[1]
        proteins[current] = ""
    else:
        proteins[current] += line

# ------------------------------------------------------------------------------------------------------
# create a dataframe of Protein ID | Sequence | EC Number

frame = pd.DataFrame(columns=["Protein ID", "Sequence", "EC Number"])
i = 0
for key, value in proteins.items():
    try:
        frame = frame.append(pd.DataFrame(data=pd.Series([[key, value, list(ECs.keys())[list(ECs.values()).index(key)]]]), columns=["Protein ID", "Sequence", "EC Number"]))
    except ValueError:
        frame = frame.append(pd.DataFrame(data=[[key, value, "NONE"]], columns=["Protein ID", "Sequence", "EC Number"]))
    i += 1
    if i % 100 == 0:
        print("Creating row %d of %d" % (i, len(proteins)))

frame.to_csv(Path.home() + "/Annotation/Data/frame", sep="\t")
