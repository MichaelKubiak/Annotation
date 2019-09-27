# ------------------------------------------------------------------------------------------------------
# imports

from scipy import sparse
from paths import DATA
import numpy as np
import prep


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

    targets = sparse.load_npz(DATA + "target_matrix.npz").todense()

    # ------------------------------------------------------------------------------------------------------
    # Print percentage of empty rows in the target matrix

    print("Percentage empty rows in target matrix before pruning: %.2f%%" % (prep.get_empty(targets)))


    prep.remove_non_family(scores, targets)
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
