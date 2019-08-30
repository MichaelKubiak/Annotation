#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to parse the hmmer results table into a matrix of scores
# ------------------------------------------------------------------------------------------------------

from paths import DATA
import re
from scipy import sparse
import numpy as np
from scipy import io


# ------------------------------------------------------------------------------------------------------
# Make empty sparse matrix

with open(DATA + "Pfam-A.hmm", encoding="utf-8") as pfamfile:
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

mat = sparse.lil_matrix((len(protein_Accessions), len(pfam_Accessions)), dtype=np.float32)


# ------------------------------------------------------------------------------------------------------
# put the searchhmm results from the table into the matrix

with open(DATA + "hmmresult_full") as infile:
    result = infile.readlines()


scores = []
for line in result:
    if not line.startswith("#"):
        splitline = list(filter(None, re.split(r"\s", line)))
        scores.append((splitline[0].split("|")[1], splitline[3], splitline[5]))

for score in scores:
    mat[protein_Accessions.index(score[0]), pfam_Accessions.index(score[1])] = score[2]

# ------------------------------------------------------------------------------------------------------
# output matrix and column/row headings
print(mat.data.nbytes)
io.mmwrite(DATA + "score_matrix", mat)
with open(DATA + "matrix_rows", "w") as rows, open(DATA + "matrix_columns", "w") as columns:
    rows.write("\n".join(protein_Accessions))
    columns.write("\n".join(pfam_Accessions))


