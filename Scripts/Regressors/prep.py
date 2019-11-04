# Module for preparing the data to be learned from
# ------------------------------------------------------------------------------------------------------
# imports

import random
import numpy as np
import re
from collections import defaultdict

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to split a dataset and target list into training and test sets

def train_test_split_sparse(data, proteins, targets, test_size, random_state):
    rearr = list(range(data.shape[0]))
    random.seed(random_state)
    random.shuffle(rearr)
    X_shuffle = data[rearr]
    X_train = X_shuffle[round(data.shape[0]*test_size):]
    X_test = X_shuffle[:round(data.shape[0]*test_size)]
    print(len(targets))
    print(max(rearr))
    y_shuffle = [targets[x] for x in rearr]
    y_train = y_shuffle[round(len(targets)*test_size):]
    y_test = y_shuffle[:round(len(targets)*test_size)]
    proteins_shuffle = [proteins[x] for x in rearr]
    proteins_train = proteins_shuffle[round(len(proteins)*test_size):]
    proteins_test = proteins_shuffle[:round(len(proteins)*test_size)]

    return X_train, X_test, y_train, y_test, proteins_train, proteins_test


# ------------------------------------------------------------------------------------------------------
# Function to return the percentage of empty rows in a matrix

def get_empty(list):
    return 100*len([list == "None"])/len(list)


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
    targets = [targets[row] for row in range(len(targets)) if rows[row]]

    return scores, proteins, pfam, targets


# ------------------------------------------------------------------------------------------------------
# Function to remove proteins with no EC number down to a specified limit

def remove_non_enzyme(scores, proteins, targets, limit):
    i = 0
    # While the ratio of empty rows in target is greater than the limit
    while (len(targets) - len([targets != "None"]))/len(targets) > limit:
        # Set a different random seed each time, so that the same numbers aren't generated repeatedly
        np.random.seed(i)
        # Generate 100,000 random integers < the number of rows
        rands = list(np.random.randint(0, len(targets), 100000, dtype=int))
        rem = []
        # if the row in targets for a number is not empty, remove it from the list
        for r in rands:
            if targets[r] != "None":
                rem.append(r)
        for x in rem:
            rands.remove(x)

        # Remove the rows from both matrices and the protein list

        # Make a mask of Trues the same length as the number of rows in scores
        mask = np.ones(scores.shape[0], dtype=bool)
        # For numbers still in the list, make that index in the mask False
        mask[rands] = False
        # Remove all rows with False from scores
        targets = [targets[pos] for pos in range(len(proteins)) if mask[pos]]
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
# Function to order ECs and convert them to numerical values - also converts the target list

def order_ECs(targets):
    ECs = []
    for protein in range(len(targets)):
        if targets[protein] != "None":
            for EC in targets[protein].split("\t"):
                if EC not in ECs:
                    ECs.append(EC)
    order = Nestdict()
    for EC in ECs:
        splitEC = EC.split(".")
        if re.search("n", splitEC[3]):
            splitEC[3] = ["n", int(splitEC[3].strip("n"))]
        else:
            splitEC[3] = ["", int(splitEC[3])]
        try:
            order[splitEC[0]][splitEC[1]][splitEC[2]].append(splitEC[3])
        except AttributeError:
            order[splitEC[0]][splitEC[1]][splitEC[2]] = [splitEC[3]]

    ECs = []
    for k in sorted(order.keys()):
        for k2 in sorted(order[k].keys()):
            for k3 in sorted(order[k][k2].keys()):
                for v4 in sorted(order[k][k2][k3]):
                    ECs.append(k + "." + k2 + "." + k3 + "." + "".join(str(x) for x in v4))

    ECs = ["None"] + ECs
    con_targets = []
    for protein in targets:
        con_targets.append([])
        for EC in protein.split("\t"):
            con_targets[-1].append(ECs.index(EC))
    return ECs, con_targets


# ------------------------------------------------------------------------------------------------------
# Nested dictionary class which creates new outer dictionaries if they are required

class Nestdict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

