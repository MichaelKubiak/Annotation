# Module for preparing the data to be learned from
# ------------------------------------------------------------------------------------------------------
# imports

import random
import numpy as np
import re


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


# ------------------------------------------------------------------------------------------------------
# Function to determine the number of proteins that have each set of ECs of a certain length

def get_repeated(EC_classes, n_repeats):
    lengths = []
    repeats = {}
    # Make the list of lengths
    for protein in EC_classes:
        if protein == "None":
            lengths.append(0)
        else:
            lengths.append(len(protein.split("\t")))
            # use that list to build a dictionary of the number of times a specific group of proteins (of length n_repeats) is present
            if lengths[-1] == n_repeats:
                if protein in repeats.keys():
                    repeats[protein] += 1
                else:
                    repeats[protein] = 1
    return lengths, repeats


# ------------------------------------------------------------------------------------------------------
# Function to make a dictionary of EC numbers with which proteins have that annotation

def get_EC_dict(enzyme):
    EC_dict = {}

    for line in enzyme:
        if line.startswith("ID"):
            current = line.split()[1]
            EC_dict[current] = []
        elif line.startswith("DR"):
            for protein in re.finditer(r"\w+,", line):
                EC_dict[current].append(protein.group().split(",")[0])

    return EC_dict


# ------------------------------------------------------------------------------------------------------
# Function to remove one of each pair of EC numbers with identical membership

def remove_duplicate_ECs(EC_classes, enzyme, ECs, targets):

    lengths, repeats = get_repeated(EC_classes, 2)
    EC_dict = get_EC_dict(enzyme)
    print("EC pair", "\t|\t", "occurences", "\t|\t", "EC1 occurences", "\t|\t", "EC2 occurences")

    rem = []
    dupe_ECs = ["# kept == removed"]
    for k in repeats.keys():
        EC_1 = k.split("\t")[0]
        EC_2 = k.split("\t")[1]
        occ_1 = len(EC_dict[EC_1])
        occ_2 = len(EC_dict[EC_2])
        print(k, "\t|\t", repeats[k], "\t|\t", occ_1, "\t|\t", occ_2)
        if repeats[k] == occ_1 == occ_2:
            print(EC_1, "cut out due to identical membership to", EC_2)
            rem.append(ECs.index(EC_2))
            dupe_ECs.append("".join((EC_1, "==", EC_2)))
    dupe_rem_EC_order = [ECs[i] for i in range(len(ECs)) if i not in rem]
    indices = [range(len(ECs))[i] for i in range(len(ECs)) if i not in rem]
    dupe_rem_targets = targets[:, indices]
    return dupe_ECs, dupe_rem_EC_order, dupe_rem_targets
