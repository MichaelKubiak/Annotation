#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to produce statistics about the data
# ------------------------------------------------------------------------------------------------------
# imports

from paths import DATA
from scipy import io
import numpy as np


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    # ------------------------------------------------------------------------------------------------------
    # Calculate percentage non enzyme

    with open(DATA + "targets") as tfile:
        targets = tfile.readlines()

    targets = list(map(str.strip, targets))

    print("Non-enzyme: %.2f%%" % (targets.count("None")/len(targets)*100))

    # ------------------------------------------------------------------------------------------------------
    # Find HMMs not represented in Swissprot

    scores = io.mmread(DATA + "score_matrix")
    # transpose the matrix so that the empty columns can be found
    invscores = scores.transpose()

    with open(DATA + "matrix_columns") as pfamfile:
        pfam = pfamfile.readlines()

    pfam = list(map(str.strip, pfam))
    # make array of booleans showing whether there are any non-zero numbers in that row of the transposed matrix
    used = np.diff(invscores.tocsr().indptr) != 0

    unused = np.where(used == False)

    unusedhmm = []
    for col in unused[0]:
        unusedhmm.append(pfam[col])

    print("Unused HMMs: number = %d %s" % (len(unused), str(unusedhmm)))

    # ------------------------------------------------------------------------------------------------------
    # Find Proteins not represented by Pfam HMMs

    with open(DATA + "matrix_rows") as spfile:
        swissprot = spfile.readlines()

    swissprot = list(map(str.strip, swissprot))

    modelled = np.diff(scores.tocsr().indptr) != 0

    unmodelled = np.where(modelled == False)

    unmodelledprotein = []
    for row in unmodelled[0]:
        unmodelledprotein.append(swissprot[row])

    print("Unrepresented proteins: number = %d %s" % (len(unmodelledprotein), str(unmodelledprotein)))

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
