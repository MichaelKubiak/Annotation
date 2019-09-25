#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# Script to teach the network to determine EC numbers
# ------------------------------------------------------------------------------------------------------
# imports

from scipy import sparse
import numpy as np
from paths import DATA
from sklearn.model_selection import train_test_split

# ------------------------------------------------------------------------------------------------------
# read files

scores = sparse.load_npz(DATA + "score_matrix.npz")
# proteins are on rows in both cases
with open(DATA + "matrix_rows") as protein_accessions, open(DATA + "matrix_columns") as pfam_accessions, open(DATA + "EC_order") as EC_order:
    proteins, pfam, ECs = protein_accessions.readlines(), pfam_accessions.readlines(), EC_order.readlines()

targets = sparse.load_npz(DATA + "target_matrix.npz").todense()


print("Percentage empty rows in target matrix before pruning: %.2f%%" % (100*np.sum(~targets.any(1))/targets.shape[0]))
# find rows with only 0s
modelled = np.diff(scores.indptr) != 0

unmodelled = np.where(modelled == False)

# remove rows and columns with only 0s
scores = scores[scores.getnnz(1) > 0][:, scores.getnnz(0) > 0]

# remove those rows from targets
targets = np.delete(targets, unmodelled[0], axis=0)

print("Percentage empty rows in target matrix after pruning: %.2f%%" % (100*np.sum(~targets.any(1))/targets.shape[0]))

proteins = list(map(str.strip, proteins))
pfam = list(map(str.strip, pfam))
ECs = list(map(str.strip, ECs))

i = 0
while np.sum(~targets.any(1))/targets.shape[0] > 0.2:
    np.random.seed(i)
    rands = list(np.random.randint(0, targets.shape[0], 100000, dtype=int))
    rem = []
    for r in rands:
        if targets[r].any():
          rem.append(r)
    for x in rem:
        rands.remove(x)

    targets = np.delete(targets, rands, axis=0)
    mask = np.ones(scores.shape[0], dtype=bool)
    mask[rands] = False
    scores = scores[mask]
    i += 1
    print("Pass %d: %.2f%%" % (i, 100*np.sum(~targets.any(1))/targets.shape[0]))
print(scores.shape)
print(targets.shape)

print("Percentage empty rows in target matrix after further pruning: %.2f%%" % (100*np.sum(~targets.any(1))/targets.shape[0]))

#X_train, X_test, y_train, y_test = train_test_split(scores, targets, test_size=.3, random_state=1, stratify=targets)

# RAM required, smaller sample test set
