#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to produce statistics about the data
# ------------------------------------------------------------------------------------------------------
# imports

from paths import DATA
from scipy import sparse
import numpy as np


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    # ------------------------------------------------------------------------------------------------------
    # Read file for percentage calculation

    with open(DATA + "targets") as tfile:
        targets = tfile.readlines()

    # ------------------------------------------------------------------------------------------------------
    # Calculate percentage non enzyme

    targets = list(map(str.strip, targets))
    print(targets.count("None"))
    print("Non-enzyme: %.2f%%" % (targets.count("None")/len(targets)*100))

    # ------------------------------------------------------------------------------------------------------
    # Read in files for finding unrepresented HMMs and proteins

    # load .npz file
    scores = sparse.load_npz(DATA + "score_matrix.npz")

    # transpose the matrix so that the empty columns can be found
    scores_t = scores.transpose()

    with open(DATA + "matrix_columns") as pfamfile, open(DATA + "matrix_rows") as spfile:
        pfam, swissprot = pfamfile.readlines(), spfile.readlines()

    # strip newline characters from the ends of list members
    pfam, swissprot = list(map(str.strip, pfam)), list(map(str.strip, pfam))

    # ------------------------------------------------------------------------------------------------------
    # Find HMMs not represented in Swissprot

    # Make array of booleans showing whether there are any non-zero numbers in that row of the transposed matrix
    used = np.diff(scores_t.tocsr().indptr) != 0

    # Unused is a list of indices at which used is False
    unused = np.where(used == False)

    # Make a list of family names from the list of indices
    unusedhmm = []
    for col in unused[0]:
        unusedhmm.append(pfam[col])

    # Print the number of members of the list followed by the list
    print("Unused HMMs: number = %d %s" % (len(unusedhmm), str(unusedhmm)))

    # ------------------------------------------------------------------------------------------------------
    # Find Proteins not represented by Pfam HMMs

    # Repeat the process above for the proteins instead of the families
    modelled = np.diff(scores.tocsr().indptr) != 0

    unmodelled = np.where(modelled == False)

    unmodelledprotein = []
    for row in unmodelled[0]:
        unmodelledprotein.append(swissprot[row])

    print("Unrepresented proteins: number = %d %s" % (len(unmodelledprotein), str(unmodelledprotein)))


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
