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


# ------------------------------------------------------------------------------------------------------
# Function to remove proteins with no EC number down to a specified limit

def remove_non_enzyme(scores, targets, limit):
    i = 0
    # While the ratio of empty rows in target is greater than the limit
    while np.sum(~targets.any(1))/targets.shape[0] > limit:
        # Set a different random seed each time, so that the same numbers aren't generated repeatedly
        np.random.seed(i)
        # Generate 100,000 random integers < the number of rows
        rands = list(np.random.randint(0, targets.shape[0], 100000, dtype=int))
        rem = []
        # if the row in targets for a number is not empty, remove it from the list
        for r in rands:
            if targets[r].any():
                rem.append(r)
        for x in rem:
            rands.remove(x)

        # Remove all rows for numbers still in the list from targets (targets is a dense matrix)
        targets = np.delete(targets, rands, axis=0)

        # Remove the same rows from scores (scores is a sparse matrix)

        # Make a mask of Trues the same length as the number of rows in scores
        mask = np.ones(scores.shape[0], dtype=bool)
        # For numbers still in the list, make that index in the mask False
        mask[rands] = False
        # Remove all rows with False from scores
        scores = scores[mask]
        i += 1
        # Tell the user how far it has got, and what the current percentage is
        print("Pass %d: %.2f%%" % (i, get_empty(targets)))
    return scores, targets


