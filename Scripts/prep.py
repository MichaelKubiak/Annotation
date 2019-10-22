# Module for preparing the data to be learned from
# ------------------------------------------------------------------------------------------------------
# imports

import random
import numpy as np


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to split a dataset and target matrix into training and test sets while BOTH ARE SPARSE

def train_test_split_sparse(data, proteins, targets, test_size, random_state):
    rearr = list(range(data.shape[0]))
    random.seed(random_state)
    random.shuffle(rearr)
    X_shuffle = data[rearr]
    X_train = X_shuffle[round(data.shape[0]*test_size):]
    X_test = X_shuffle[:round(data.shape[0]*test_size)]
    y_shuffle = targets[rearr]
    y_train = y_shuffle[round(targets.shape[0]*test_size):]
    y_test = y_shuffle[:round(targets.shape[0]*test_size)]
    proteins_shuffle = [proteins[x] for x in rearr]
    proteins_train = proteins_shuffle[round(len(proteins)*test_size):]
    proteins_test = proteins_shuffle[:round(len(proteins)*test_size)]

    return X_train, X_test, y_train, y_test, proteins_train, proteins_test


# ------------------------------------------------------------------------------------------------------
# Function to return the percentage of empty rows in a matrix

def get_empty(matrix):
    return 100*np.sum(~matrix.todense().any(1))/matrix.shape[0]


# ------------------------------------------------------------------------------------------------------
# Function to remove proteins that are not hit by any pfam hmms from both matrices

def remove_non_family(scores, proteins, pfam, targets):

    # find rows and columns with non-zero values
    rows = scores.getnnz(1) > 0
    columns = scores.getnnz(0) > 0

    # slice the score matrix to only contain those rows and columns
    scores = scores[rows][:, columns]
    proteins = [proteins[row] for row in range(len(proteins)) if rows[row]]
    pfam = [pfam[column] for column in range(len(pfam)) if columns[column]]

    # remove those rows from targets
    targets = targets[rows]

    return scores, proteins, pfam, targets


# ------------------------------------------------------------------------------------------------------
# Function to remove proteins with no EC number down to a specified limit

def remove_non_enzyme(scores, proteins, targets, limit):
    i = 0
    # While the ratio of empty rows in target is greater than the limit
    while (targets.shape[0] - np.sum(targets.getnnz(1)))/targets.shape[0] > limit:
        # Set a different random seed each time, so that the same numbers aren't generated repeatedly
        np.random.seed(i)
        # Generate 100,000 random integers < the number of rows
        rands = list(np.random.randint(0, targets.shape[0], 100000, dtype=int))
        rem = []
        # if the row in targets for a number is not empty, remove it from the list
        for r in rands:
            if targets[r].getnnz():
                rem.append(r)
        for x in rem:
            rands.remove(x)

        # Remove the rows from both matrices and the protein list

        # Make a mask of Trues the same length as the number of rows in scores
        mask = np.ones(scores.shape[0], dtype=bool)
        # For numbers still in the list, make that index in the mask False
        mask[rands] = False
        # Remove all rows with False from scores
        targets = targets[mask]
        scores = scores[mask]
        proteins = [proteins[pos] for pos in range(len(proteins)) if mask[pos]]
        i += 1
        # Tell the user how far it has got, and what the current percentage is
        print("Pass %d: %.2f%%" % (i, get_empty(targets)))
    return scores, proteins, targets


