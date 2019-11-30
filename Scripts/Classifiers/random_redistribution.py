#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to randomly redistribute the rows of the expected output matrix as a baseline for accuracy
# ------------------------------------------------------------------------------------------------------
# Imports

from Classifiers import prep,train_model as tm,test_harness as th
from Classifiers.test_data import create_test
import random
import numpy as np

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Test the model

def test_rearrange(random_state, y_test, ECs):

    rearr = list(range(y_test.shape[0]))
    random.seed(random_state)
    random.shuffle(rearr)
    pred = y_test[rearr]
    return th.get_metrics(pred.todense(), y_test.todense(), ECs)


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    args = tm.arguments("random redistribution", "random.clf")

    # ------------------------------------------------------------------------------------------------------
    # Read files

    proteins, pfam, ECs, scores, targets = tm.read_files(args)

    # ------------------------------------------------------------------------------------------------------
    # make a test dataset

    scores, proteins, pfam, targets, ECs = create_test(scores, proteins, pfam, targets, ECs)

    # ------------------------------------------------------------------------------------------------------
    # Remove proteins with no pfam hits - nothing happens with test set

    print("Percentage empty rows in target matrix before pruning: %.2f%%"%(prep.get_empty(targets)))

    scores, proteins, pfam, targets = prep.remove_non_family(scores,proteins,pfam,targets)

    print("Percentage empty rows in target matrix after pruning: %.2f%%"%(prep.get_empty(targets)))

    # ------------------------------------------------------------------------------------------------------
    # Remove non-enzyme proteins down to a limit

    limit = 0.2

    scores, proteins, targets = prep.remove_non_enzyme(scores, proteins, targets, limit)

    print("Percentage empty rows in target matrix after removal of empty rows down to %d: %.2f%%" % (limit, prep.get_empty(targets)))

    #test method
    test_scores = []
    for i in range(10):
        print("rand =", i)
        # Make learning and test datasets
        X_train, X_test, y_train, y_test, proteins_train, proteins_test = prep.train_test_split_sparse(scores,proteins,targets,test_size=0.3,random_state=i)
        test_scores.append(test_rearrange(i, y_test, ECs))
    test_scores = np.array(test_scores)
    print("Total mean accuracy:", np.mean(test_scores[:, 0]))
    print("Total mean sensitivity:", np.mean(test_scores[:, 1]))
    print("Total mean specificity:", np.mean(test_scores[:, 2]))
    print("Total mean precision:", np.mean(test_scores[:, 3]))
    # print("Total mean F1 score:", (2*mean(test_scores[:, 1])*mean(test_scores[:, 3])/(mean(test_scores[:, 1] + mean(test_scores[:, 3])))))
    print(test_scores[0])

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
