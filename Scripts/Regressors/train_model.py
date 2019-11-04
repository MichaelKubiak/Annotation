# ------------------------------------------------------------------------------------------------------
# A module containing methods to build and train a classifier
# ------------------------------------------------------------------------------------------------------
# Imports

import argparse
from paths import DATA
from scipy import sparse


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Set up argparse

def arguments(descriptor, output):
    parser = argparse.ArgumentParser(description="A script to generate a %s to determine the EC numbers of proteins" % descriptor)
    parser.add_argument("-p", "--path", default=DATA, help="Path to the folder to be used for i/o")
    parser.add_argument("-s", "--score", default="score_matrix.npz", help="File name of the score matrix")
    parser.add_argument("-t", "--targets", default="targets", help="File name of the target list")
    parser.add_argument("-a", "--protein_accessions", default="matrix_rows", help="File name of the list of protein accessions")
    parser.add_argument("-f", "--pfam_accessions", default="matrix_columns", help="File name of the list of pfam accessions")
    parser.add_argument("-o", "--output", default=output, help="File name for output of the model")
    return parser.parse_args()


# ------------------------------------------------------------------------------------------------------
# Read in necessary files

def read_files(args):
    scores = sparse.load_npz(args.path+args.score)
    # Proteins are on rows in both cases
    with open(args.path+args.protein_accessions) as protein_accessions, open(args.path+args.pfam_accessions) as pfam_accessions,\
            open(args.path+args.targets) as ECs:
        proteins, pfam, targets = protein_accessions.readlines(), pfam_accessions.readlines(), ECs.readlines()
    proteins = list(map(str.strip, proteins))
    pfam = list(map(str.strip, pfam))
    targets = list(map(str.strip, targets))
    return proteins, pfam, scores, targets
