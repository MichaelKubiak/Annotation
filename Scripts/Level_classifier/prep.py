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

    shuffle = rearrange(current, random_state)

    train = shuffle[round(len(current)*test_size):]
    test = shuffle[:round(len(current)*test_size)]

    unique_targets = set([target.split(".")[level] for group in [target_set.split("\t") for target_set in [targets[x] for x in current]] for target in group])
    return train, test, unique_targets


# ------------------------------------------------------------------------------------------------------
# Function to rearrange the data and targets

def rearrange(current, random_state):

    rearr = list(range(len(current)))
    random.seed(random_state)
    random.shuffle(rearr)
    return [current[x] for x in rearr]


# ------------------------------------------------------------------------------------------------------
# Function to get the current level of each classification

def get_level(targets, level):

    return ["\t".join([num.split(".")[level] for num in nums]) for nums in [target.split("\t") for target in targets]]
