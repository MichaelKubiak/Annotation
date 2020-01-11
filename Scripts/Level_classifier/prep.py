# Module for preparing the data to be learned from
# ------------------------------------------------------------------------------------------------------
# imports

import random

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to split a dataset and target list into training and test sets

def train_test_split_sparse(data, targets, test_size, level=1, prev_class="None", random_state=0):
    if prev_class is not "None":
        current = [x for x in len(targets) if targets[x].split(".")[:level] == prev_class]
    else:
        current = range(len(targets))
