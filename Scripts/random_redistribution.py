#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to randomly redistribute the rows of the expected output matrix as a baseline for accuracy
# ------------------------------------------------------------------------------------------------------
# Imports

import train_model as tm
import prep
from test_data import create_test
from sklearn.model_selection import train_test_split
import test_harness as th
import random
import numpy as np
from statistics import mean, pstdev


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Test the model

def test_rearrange(random_state, y_test):

    rearr = list(range(y_test.shape[0]))
    random.seed(random_state)
    random.shuffle(rearr)
    pred = y_test[rearr]
    return th.get_metrics(pred, y_test)


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

    #test method
    test_scores = []
    for i in range(10):
        print("rand =", i)
        # Make learning and test datasets
        X_train, X_test, y_train, y_test = train_test_split(scores, targets, test_size=0.3, random_state=i)
        test_scores.append(test_rearrange(i, y_test))
    test_scores = np.array(test_scores)
    print("Total mean accuracy:", mean(test_scores[:, 0]))
    print("Total mean sensitivity:", mean(test_scores[:, 1]))
    print("Total mean specificity:", mean(test_scores[:, 2]))
    print("Total mean precision:", mean(test_scores[:, 3]))
    print("Total mean F1 score:", (2*mean(test_scores[:, 1])*mean(test_scores[:, 3])/(mean(test_scores[:, 1] + mean(test_scores[:, 3])))))


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
