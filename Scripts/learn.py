# ------------------------------------------------------------------------------------------------------
# Script to teach the network to determine EC numbers
# ------------------------------------------------------------------------------------------------------
# imports

from scipy import io
import numpy as np
from paths import DATA
from sklearn.model_selection import train_test_split
# ------------------------------------------------------------------------------------------------------
# read files

scores = io.mmread(DATA + "score_matrix").tocsr()
# proteins are on rows in both cases
with open(DATA + "matrix_rows") as protein_accessions, open(DATA + "matrix_columns") as pfam_accessions, open(DATA + "EC_order") as EC_order:
    proteins, pfam, ECs = protein_accessions.readlines(), pfam_accessions.readlines(), EC_order.readlines()

targets = io.mmread(DATA + "target_matrix.mtx").astype(np.bool).todense()

# find rows with only 0s
modelled = np.diff(scores.tocsr().indptr) != 0

unmodelled = np.where(modelled == False)

# remove rows and columns with only 0s
scores = scores[scores.getnnz(1) > 0][:, scores.getnnz(0) > 0]

# remove those rows from targets
targets = np.delete(targets, unmodelled[0], axis=0)

print(scores.shape)
print(targets.shape)

proteins = list(map(str.strip, proteins))
pfam = list(map(str.strip, pfam))
ECs = list(map(str.strip, ECs))

X_train, X_test, y_train, y_test = train_test_split(scores, targets, test_size=.3, random_state=1, stratify=targets)

# TODO: Needs running on Spectre (RAM), write whole thing before testing
