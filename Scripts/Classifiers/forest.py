#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to build a random forest classifier
# ------------------------------------------------------------------------------------------------------
# Imports

from Classifiers import prep,train_model as tm,test_harness as th
from Classifiers.test_data import create_test
from sklearn.ensemble import RandomForestClassifier
import joblib
import numpy as np
from statistics import mean

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Train a random forest classifier

def train_forest(scores, proteins, targets, i):

    # Make learning and test datasets
    X_train, X_test, y_train, y_test, proteins_train, proteins_test = prep.train_test_split_sparse(scores,proteins,targets,test_size=0.3,random_state=i)
    # Make the classifier
    model = RandomForestClassifier(random_state=i, n_estimators=1000)
    # Train the classifier on the data
    model.fit(X_train, y_train.todense())
    return X_test, model, y_test.todense(), proteins_test

# ------------------------------------------------------------------------------------------------------

def main():

    args = tm.arguments("random forest", "forest.clf")

    # ------------------------------------------------------------------------------------------------------
    # Read files

    proteins, pfam, ECs, scores, targets = tm.read_files(args)

    # ------------------------------------------------------------------------------------------------------
    # make a test dataset
    print(scores.shape)
    scores, proteins, pfam, targets, ECs = create_test(scores, proteins, pfam, targets, ECs)
    print(scores.shape)
    # ------------------------------------------------------------------------------------------------------
    # Remove proteins with no pfam hits - nothing happens with test set

    print("Percentage empty rows in target matrix before pruning: %.2f%%"%(prep.get_empty(targets)))

    scores, proteins, pfam, targets = prep.remove_non_family(scores,proteins,pfam,targets)

    print("Percentage empty rows in target matrix after pruning: %.2f%%"%(prep.get_empty(targets)))

    # ------------------------------------------------------------------------------------------------------
    # Remove non-enzyme proteins down to a limit

    limit = 0.2

    scores, proteins, targets = prep.remove_non_enzyme(scores, proteins, targets, limit)

    print("Percentage empty rows in target matrix after removal of empty rows down to %.2f: %.2f%%" % (limit, prep.get_empty(targets)))
    # ------------------------------------------------------------------------------------------------------
    # test method
    test_scores = []
    for i in range(10):
        print("rand =", i)
        X_test, forest, y_test, proteins_test = train_forest(scores, proteins, targets, i)
        test_scores.append(th.test_model(X_test, forest, y_test, ECs))
    test_scores = np.array(test_scores)
    print("Total mean accuracy:", np.mean(test_scores[:, 0]))
    print("Total mean sensitivity:", np.mean(test_scores[:, 1]))
    print("Total mean specificity:", np.mean(test_scores[:, 2]))
    print("Total mean precision:", np.mean(test_scores[:, 3]))
    print("Total mean F1 score:", (2*mean(test_scores[:, 1])*mean(test_scores[:, 3])/(mean(test_scores[:, 1] + mean(test_scores[:, 3])))))
    # Output the classifier as a pickle using joblib
    # X_test, forest, y_test, proteins_test = train_forest(scores, proteins, targets, 1)
    # joblib.dump(forest, args.path + args.output)


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
