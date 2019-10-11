#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to build a random forest classifier
# ------------------------------------------------------------------------------------------------------
# Imports

from sklearn.model_selection import train_test_split
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
# Train a random forest classifier

def train_forest(scores, targets, i, classifier):

    # Make learning and test datasets
    X_train, X_test, y_train, y_test = train_test_split(scores, targets, test_size=0.3, random_state=i) # use skmultilearn?
    # Make the classifier
    model = classifier(random_state=i, n_estimators=1000, n_jobs=-1)
    # Train the classifier on the data
    model.fit(X_train, y_train)
    return X_test, model, y_test

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
    # test method
    test_scores = []
    for i in range(10):
        print("rand =", i)
        X_test, forest, y_test = train_forest(scores, targets, i)
        test_scores.append(th.test_model(X_test, forest, y_test))
    test_scores = np.array(test_scores)
    print("Total mean accuracy:", mean(test_scores[:, 0]))
    print("Total mean sensitivity:", mean(test_scores[:, 1]))
    print("Total mean specificity:", mean(test_scores[:, 2]))
    print("Total mean precision:", mean(test_scores[:, 3]))
    print("Total mean F1 score:", (2*mean(test_scores[:, 1])*mean(test_scores[:, 3])/(mean(test_scores[:, 1] + mean(test_scores[:, 3])))))
    # Output the classifier as a pickle using joblib
    X_test, forest, y_test = train_forest(scores, targets, 1, RandomForestClassifier)
    joblib.dump(forest, args.path + args.forest_output)



# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
