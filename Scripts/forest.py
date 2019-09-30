#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to build a random forest classifier
# ------------------------------------------------------------------------------------------------------
# Imports

from scipy import sparse
from paths import DATA
import numpy as np
import prep
from test_data import create_test
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import joblib


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

    targets = sparse.load_npz(DATA + "target_matrix.npz")#.todense() add when not using test dataset

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
    # make learning and test datasets

    X_train, X_test, y_train, y_test = train_test_split(scores, targets, test_size=0.3, random_state=1, stratify=targets)
    forest = RandomForestClassifier(random_state=1, n_estimators=10)

    forest.fit(X_train, y_train)
    Z = forest.predict(X_test)
    equality = (Z == y_test)
    print("Accuracy of model: %.2f%%" % (100*np.count_nonzero(equality)/y_test.size))
    joblib.dump(forest, DATA + "forest.clf")


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
