#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to build a random forest classifier
# ------------------------------------------------------------------------------------------------------
# Imports

import argparse
from scipy import sparse
from paths import DATA
import prep
from test_data import create_test
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import joblib
import test_harness as th
import numpy as np


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Train a random forest classifier
def train_forest(scores, targets, i):

    # Make learning and test datasets
    X_train, X_test, y_train, y_test = train_test_split(scores, targets, test_size=0.3, random_state=i) # use skmultilearn?
    # Make the classifier
    forest = RandomForestClassifier(random_state=i, n_estimators=1000, n_jobs=-1)
    # Train the classifier on the data
    forest.fit(X_train, y_train)
    return X_test, forest, y_test


# ------------------------------------------------------------------------------------------------------
# Test the model

def test_forest(X_test, forest, y_test):

    pred = th.predict(X_test, forest)
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

    parser = argparse.ArgumentParser(description="A script to generate a random forest to determine the EC numbers of proteins")
    parser.add_argument("-p", "--path", default=DATA, help="Path to the folder to be used for i/o")
    parser.add_argument("-s", "--score", default="score_matrix.npz", help="File name of the score matrix")
    parser.add_argument("-t", "--targets", default="target_matrix.npz", help="File name of the target matrix")
    parser.add_argument("-a", "--protein_accessions", default="matrix_rows", help="File name of the list of protein accessions")
    parser.add_argument("-f", "--pfam_accessions", default="matrix_columns", help="File name of the list of pfam accessions")
    parser.add_argument("-E", "--EC_numbers", default="EC_order", help="File name of the list of EC numbers")
    parser.add_argument("-o", "--forest_output", default="forest.clf", help="File name for output of the forest")
    args = parser.parse_args()

    # ------------------------------------------------------------------------------------------------------
    # Read files

    scores = sparse.load_npz(args.path + args.score)
    # Proteins are on rows in both cases
    with open(args.path + args.protein_accessions) as protein_accessions, open(args.path + args.pfam_accessions) as pfam_accessions, open(args.path + args.EC_numbers) as EC_order:
        proteins, pfam, ECs = protein_accessions.readlines(), pfam_accessions.readlines(), EC_order.readlines()

    proteins = list(map(str.strip, proteins))
    pfam = list(map(str.strip, pfam))
    ECs = list(map(str.strip, ECs))

    targets = sparse.load_npz(args.path + args.targets)#.todense() # add when not using test dataset

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
        X_test, forest, y_test = train_forest(scores, targets, i)
        test_scores.append(test_forest(X_test, forest, y_test))
    test_scores = np.array(test_scores)
    print("mean accuracy:", sum(test_scores[:, 0])/len(test_scores[:, 0]))
    print("mean sensitivity:", sum(test_scores[:, 1])/len(test_scores[:, 1]))
    print("mean specificity:", sum(test_scores[:, 2])/len(test_scores[:, 2]))
    # Output the classifier as a pickle using joblib
    X_test, forest, y_test = train_forest(scores, targets, 1)
    joblib.dump(forest, args.path + args.forest)


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
