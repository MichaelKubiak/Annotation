# Module for preparing the data to be learned from
# ------------------------------------------------------------------------------------------------------
# imports

import random


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to split a dataset and target list into training and test sets

def train_test_split(data, targets, test_size,  prev_class, level=0, random_state=0):
    if level != 0:
        current = [x for x in range(len(targets)) if prev_class in [".".join(target.split(".")[:level]) for target in targets[x].split("\t")]]
    else:
        current = range(len(targets))

    X_shuffle, y_shuffle = rearrange(data[current], [targets[x] for x in current], level, random_state)

    X_train = X_shuffle[round(len(current)*test_size):]
    X_test = X_shuffle[:round(len(current)*test_size)]

    y_train = y_shuffle[round(len(current)*test_size):]
    y_test = y_shuffle[:round(len(current)*test_size)]
    return X_train, X_test, y_train, y_test


# ------------------------------------------------------------------------------------------------------
# Function to rearrange the data and targets

def rearrange(data, targets, level, random_state):

    rearr = list(range(data.shape[0]))
    random.seed(random_state)
    random.shuffle(rearr)
    X_shuffle = data[rearr]
    # the levelth number from each \t separated group in the list
    y = ["\t".join([num.split(".")[level] for num in nums]) for nums in [target.split("\t") for target in targets]]

    y_shuffle = [y[i] for i in rearr]

    return X_shuffle, y_shuffle
