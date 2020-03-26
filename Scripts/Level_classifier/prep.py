# Module for preparing the data to be learned from
# ------------------------------------------------------------------------------------------------------
# imports

import random


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to split a dataset and target list into training and test sets

def train_test_split(targets, test_size,  prev_class, level=0, random_state=0):
    if level != 0:
        current = [x for x in range(len(targets)) if prev_class in [".".join(target.split(".")[:level]) for target in targets[x].split("\t")]]
    else:
        current = range(len(targets))

    pairs = rearrange(targets, current, level, random_state)

    pairs_train = pairs[round(len(current)*test_size):]
    pairs_test = pairs[:round(len(current)*test_size)]

    return pairs_train, pairs_test


# ------------------------------------------------------------------------------------------------------
# Function to rearrange the data and targets

def rearrange(targets, current, level, random_state):

    # create tuples of (matrix row number, classification)
    pairs = []
    for i in range(len(targets)) and current:
        for target in targets[i].split("\t"):
            pairs.append((i, target.split(".")[level]))

    # randomly rearrange the tuples in the list
    rearr = list(range(len(pairs)))
    random.seed(random_state)
    random.shuffle(rearr)
    pairs_shuffle = [pairs[i] for i in rearr]

    return pairs_shuffle
