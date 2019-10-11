#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to build a random forest classifier
# ------------------------------------------------------------------------------------------------------
# Imports

import argparse
from scipy import sparse
import prep
from test_data import create_test
from sklearn.ensemble import RandomForestClassifier
import joblib
import test_harness as th
import numpy as np
from statistics import mean, pstdev
import train_model as tm


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    args = tm.arguments("random forest", "forest.clf")

    # ------------------------------------------------------------------------------------------------------
    # Read files

    proteins, pfam, ECs, scores, targets = tm.read_files(args)

    # ------------------------------------------------------------------------------------------------------
    # make a test dataset

    scores, targets = create_test(scores, pfam, targets)

    # ------------------------------------------------------------------------------------------------------
    # Remove proteins with no pfam hits - nothing happens with test set

    print("Percentage empty rows in target matrix before pruning: %.2f%%" % (prep.get_empty(targets)))

    scores, targets = prep.remove_non_family(scores, targets)

    print("Percentage empty rows in target matrix after pruning: %.2f%%" % (prep.get_empty(targets)))

    # ------------------------------------------------------------------------------------------------------
    # Remove non-enzyme proteins down to a limit

    # limit = 0.2
    #
    # scores, targets = prep.remove_non_enzyme(scores, targets, limit)
    #
    # print("Percentage empty rows in target matrix after removal of empty rows down to %d: %.2f%%" % (limit, prep.get_empty(targets)))
    # ------------------------------------------------------------------------------------------------------
    # Test method

    th.test_method(scores, targets, RandomForestClassifier)
    # ------------------------------------------------------------------------------------------------------
    # Output a classifier as a pickle using joblib - To be changed later
    X_test, forest, y_test = tm.train_model(scores, targets, 1, RandomForestClassifier)
    joblib.dump(forest, args.path + args.output)


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
