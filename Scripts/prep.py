# Module for preparing the data to be learned from
# ------------------------------------------------------------------------------------------------------
# imports

import numpy as np


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to return the percentage of empty rows in a matrix

def get_empty(matrix):
    return 100*np.sum(~matrix.any(1))/matrix.shape[0]


# ------------------------------------------------------------------------------------------------------
# Function to remove proteins that are not hit by any pfam hmms from both matrices
def remove_non_family(scores, targets):
    # Find rows with only 0s
    modelled = np.diff(scores.indptr) != 0
    unmodelled = np.where(modelled == False)
    # remove rows and columns with only 0s
    scores = scores[scores.getnnz(1) > 0][:, scores.getnnz(0)>0]
    # remove those rows from targets
    targets = np.delete(targets, unmodelled[0], axis=0)
    return scores, targets


print("Percentage empty rows in target matrix after pruning: %.2f%%" % (100*np.sum(~targets.any(1))/targets.shape[0]))


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
