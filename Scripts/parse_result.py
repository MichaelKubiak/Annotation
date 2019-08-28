#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to parse the hmmer results table into a matrix of scores
# ------------------------------------------------------------------------------------------------------

from pathlib import Path
import re
import numpy as np
import pandas
# ------------------------------------------------------------------------------------------------------
# Useful variables

HOME = str(Path.home())
DATA = HOME + "/Annotation/Data/"

# ------------------------------------------------------------------------------------------------------
# Make empty matrix

with open(DATA+ "Pfam-A.hmm") as pfamfile:
    pfam = pfamfile.readlines()

pfam_Accessions = []
for line in pfam:
    if line.startswith("ACC"):
        pfam_Accessions.append(re.split(r"\s", line)[3])

with open(DATA + "uniprot_sprot.fasta") as fastafile:
    fasta = fastafile.readlines()

protein_Accessions = []
for line in fasta:
    if line.startswith(">"):
        protein_Accessions.append(line.split("|")[1])

zeros = np.zeros(shape=(len(pfam_Accessions), len(protein_Accessions)))
frame = pandas.DataFrame(zeros, columns=pfam_Accessions, index=protein_Accessions)

print(frame)

# ------------------------------------------------------------------------------------------------------
with open(DATA + "hmmresult_56000vs1800") as infile:
    result = infile.readlines()

for line in result:
    if not line.startswith("#"):
        splitline = re.split(r"\s", line)



