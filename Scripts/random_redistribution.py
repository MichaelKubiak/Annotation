#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to randomly redistribute the rows of the expected output matrix as a baseline for accuracy
# ------------------------------------------------------------------------------------------------------
# Imports

from scipy import sparse
from paths import DATA
import prep
from test_data import create_test
from sklearn.model_selection import train_test_split
import test_harness as th
import random
import numpy as np


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Test the model

def test_rearrange(y_test):

    rearr = list(range(y_test.shape[0]))
    random.shuffle(rearr)
    pred = y_test[rearr]
    score = th.get_accuracy(pred, y_test)
    print("Accuracy of model: %.2f%%" % score)
    positives, negatives = th.get_numbers(pred, y_test)
    print("%d were false positives and %d were false negatives" % (positives[1], negatives[1]))
    sensitivity = 100*positives[0]/(positives[0] + negatives[1])
    print("Sensitivity of model: %.2f%%" % sensitivity)
    specificity = 100*negatives[0]/(negatives[0] + positives[1])
    print("Specificity of model: %.2f%%" % specificity)
    return score, sensitivity, specificity


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    # ------------------------------------------------------------------------------------------------------
    # Read files

    scores = sparse.load_npz(DATA + "score_matrix.npz")
    # Proteins are on rows in both cases
    with open(DATA + "matrix_rows") as protein_accessions, open(DATA + "matrix_columns") as pfam_accessions, open(DATA + "EC_order") as EC_order:
        proteins, pfam, ECs = protein_accessions.readlines(), pfam_accessions.readlines(), EC_order.readlines()

    proteins = list(map(str.strip, proteins))
    pfam = list(map(str.strip, pfam))
    ECs = list(map(str.strip, ECs))

    targets = sparse.load_npz(DATA + "target_matrix.npz")#.todense() # add when not using test dataset

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
        test_scores.append(test_rearrange(y_test))
    test_scores = np.array(test_scores)
    print("mean accuracy:", sum(test_scores[:,0])/len(test_scores[:,0]))
    print("mean sensitivity:", sum(test_scores[:,1])/len(test_scores[:,1]))
    print("mean specificity:", sum(test_scores[:,2])/len(test_scores[:,2]))


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
