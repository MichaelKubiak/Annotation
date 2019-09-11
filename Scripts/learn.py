# ------------------------------------------------------------------------------------------------------
# Script to teach the network to determine EC numbers
# ------------------------------------------------------------------------------------------------------
# imports

from scipy import io
from paths import DATA
from sklearn.model_selection import train_test_split
# ------------------------------------------------------------------------------------------------------
# read files

scores = io.mmread(DATA + "score_matrix")
# proteins are on rows in both cases
with open(DATA + "matrix_rows") as protein_accessions, open(DATA + "matrix_columns") as pfam_accessions, open(DATA + "EC_order") as EC_order:
    proteins, pfam, ECs = protein_accessions.readlines(), pfam_accessions.readlines(), EC_order.readlines()

targets = io.mmread(DATA + "target_matrix")

proteins = list(map(str.strip, proteins))
pfam = list(map(str.strip, pfam))
ECs = list(map(str.strip, ECs))

print(len(proteins))
print(scores.shape)
print(targets.shape)




