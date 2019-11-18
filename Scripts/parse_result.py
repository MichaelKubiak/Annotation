#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to parse the hmmer results table into a matrix of scores
# ------------------------------------------------------------------------------------------------------
#Imports

import argparse
import re
from scipy import sparse
import numpy as np


# ------------------------------------------------------------------------------------------------------
def path_arg():
    from paths import DATA
    parser = argparse.ArgumentParser(description="One of a set of scripts to parse the result of a searchhmm to be used in the model")
    parser.add_argument("-p", "--path", default=DATA, help="Path to the folder to be used for i/o")
    return parser.parse_args()


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    DATA = path_arg().path
    # ------------------------------------------------------------------------------------------------------
    # Read files

    with open(DATA + "Pfam-A.hmm", encoding="utf-8") as pfamfile:
        pfam = pfamfile.readlines()

    with open(DATA + "uniprot_sprot.fasta") as fastafile:
        fasta = fastafile.readlines()

    with open(DATA + "hmmresult_full") as infile:
        result = infile.readlines()

    # ------------------------------------------------------------------------------------------------------
    # Make empty sparse matrix

    # Make a list of all Pfam accessions from the HMM file
    pfam_Accessions = []
    for line in pfam:
        if line.startswith("ACC"):
            pfam_Accessions.append(re.split(r"\s", line)[3])

    # Make a list of all protein accessions in swissprot from the fasta file
    protein_Accessions = []
    for line in fasta:
        if line.startswith(">"):
            protein_Accessions.append(line.split("|")[1])

    # Create an empty, sparse matrix with the proportions of the lists that were made
    mat = sparse.lil_matrix((len(protein_Accessions), len(pfam_Accessions)), dtype=np.float32)

    # ------------------------------------------------------------------------------------------------------
    # Put the searchhmm results from the table into the matrix

    # Create a list of tuples that contains the accession number of a protein, the accession number of a pfam family, and the hmmer score
    scores = []
    for line in result:
        if not line.startswith("#"):
            splitline = list(filter(None, re.split(r"\s", line)))
            scores.append((splitline[0].split("|")[1], splitline[3], splitline[5]))

    # Put the score from each tuple into the matrix cell of the appropriate protein and HMM
    for score in scores:
        mat[protein_Accessions.index(score[0]), pfam_Accessions.index(score[1])] = score[2]

    # ------------------------------------------------------------------------------------------------------
    # output matrix and column/row headings

    # save the sparse matrix as a .npz file
    sparse.save_npz(DATA + "score_matrix", mat.tocsr())

    with open(DATA + "matrix_rows", "w") as rows, open(DATA + "matrix_columns", "w") as columns:
        rows.write("\n".join(protein_Accessions))
        columns.write("\n".join(pfam_Accessions))


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()

