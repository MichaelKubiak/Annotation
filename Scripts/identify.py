#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to identify Swissprot entries with EC numbers
# ------------------------------------------------------------------------------------------------------
# Imports

from parse_result import path_arg
import re
import numpy as np
from scipy import sparse


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to make a dictionary of protein IDs and the EC numbers of those proteins

def get_ECs(enzyme):
    ECs = {}
    current = ""
    for line in enzyme:
        if line.startswith("ID"):
            current = line.split()[1]
        elif line.startswith("DR"):
            for protein in re.finditer(r"\w+,", line):
                if protein.group().split(",")[0] in ECs:
                    ECs[protein.group().split(",")[0]].append(current)
                else:
                    ECs[protein.group().split(",")[0]] = [current]
    return ECs


# ------------------------------------------------------------------------------------------------------
# Function to get a list of the EC numbers for each protein in the same order as they are in the score matrix

def get_targetlist(ECs, order):
    targetlist = []
    for accession in order:
        accession = accession.strip("\n")
        try:
            targetlist.append(ECs[accession])
        except KeyError:
            targetlist.append("None")
    return targetlist


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    DATA = path_arg().path
    # ------------------------------------------------------------------------------------------------------
    # Read in files

    with open(DATA + "enzyme.dat") as efile:
        enzyme = efile.readlines()

    with open(DATA+"matrix_rows") as rows:
        order = rows.readlines()

    # ------------------------------------------------------------------------------------------------------
    # Get a dictionary of EC numbers and the IDs of proteins of that classification

    ECs = get_ECs(enzyme)

    # ------------------------------------------------------------------------------------------------------
    # Order list of EC numbers from row file

    targetlist = get_targetlist(ECs, order)

    # ------------------------------------------------------------------------------------------------------
    # Output a list of EC numbers which matches the order of the proteins in the matrix

    outstring = ""
    for target in targetlist:
        # separate the EC lists of proteins by a newline (splits sometimes leave extra empty strings when splitting on whitespace characters)
        if outstring != "":
            outstring += "\n"
        # add the EC list for the current protein, joining multiples by a tab character
        if target == "None":
            outstring += target
        else:
            outstring += "\t".join(target)

    with open(DATA + "targets", "w") as targets:
        targets.write(outstring)

    # ------------------------------------------------------------------------------------------------------
    # Output a sparse matrix with proteins as rows, and EC numbers as columns

    # Split the list into targets for each protein
    targets = outstring.split("\n")

    # Remove all repeated targets
    uniquetargets = np.unique(targets)

    # Find targets with multiple constituents
    rem = []
    for target in uniquetargets:
        if re.search("\t", target):
            rem.append(target)
            splittarget = target.split("\t")
            # Add the constituents to the list if they are not already present
            for sep in splittarget:
                if not sep in uniquetargets:
                    uniquetargets = np.append(uniquetargets, sep)

    # Remove the multiple constituent targets that have been found
    for r in rem:
        uniquetargets = uniquetargets[uniquetargets != r]
    # Remove useless targets from the list from the list of targets
    uniquetargets = uniquetargets[uniquetargets != "None"]
    uniquetargets = uniquetargets[uniquetargets != ""]
    uniquetargets = uniquetargets.tolist()

    # Create an empty, sparse boolean matrix of the correct dimensions using the scipy.sparse module
    targetmatrix = sparse.lil_matrix((len(targets), len(uniquetargets)), dtype=np.bool)

    # Fill in the sparse matrix with targets for each protein
    for i in range(len(targets)):
        splittarget = targets[i].split("\t")
        for s in splittarget:
            if s != "None":
                # Add true at the correct positions of the matrix
                targetmatrix[i, uniquetargets.index(s)] = True

    # Output the matrix as a .npz file, along with the list of targets that are in the correct order
    sparse.save_npz(DATA + "target_matrix", targetmatrix.tocsr())

    with open(DATA + "EC_order", "w") as ECfile:
        for target in uniquetargets:
            ECfile.write(target + "\n")


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
