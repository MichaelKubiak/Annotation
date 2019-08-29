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

