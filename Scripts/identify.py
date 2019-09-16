#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to identify Swissprot entries with EC numbers
# ------------------------------------------------------------------------------------------------------

from paths import DATA
import re
import numpy as np
from scipy import sparse


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# function to make a dictionary of EC numbers and the IDs of proteins of that classification

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
# function to get a list of the EC numbers for each protein in the same order as they are in the score matrix

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

    # ------------------------------------------------------------------------------------------------------
    # read in files

    with open(DATA + "enzyme.dat") as efile:
        enzyme = efile.readlines()

    # ------------------------------------------------------------------------------------------------------
    # get a dictionary of EC numbers and the IDs of proteins of that classification

    ECs = get_ECs(enzyme)

    # ------------------------------------------------------------------------------------------------------
    # read in row file from parse_result.py to order list of EC numbers

    with open(DATA + "matrix_rows") as rows:
        order = rows.readlines()

    targetlist = get_targetlist(ECs, order)

    # ------------------------------------------------------------------------------------------------------
    # output a list of EC numbers which matches the order of the proteins in the matrix

    outstring = ""
    for target in targetlist:
        if outstring != "":
            outstring += "\n"
        if target == "None":
            outstring += target
        else:
            outstring += "\t".join(target)

    with open(DATA + "targets", "w") as targets:
        targets.write(outstring)

    # ------------------------------------------------------------------------------------------------------
    # output a sparse matrix with proteins as rows, and EC numbers as columns
    # maybe unnecessary
    targets = outstring.split("\n")
    uniquetargets = np.unique(targets)

    rem = []
    for target in uniquetargets:
        if re.search("\t", target):
            rem.append(target)
            splittarget = target.split("\t")
            for sep in splittarget:
                if not sep in uniquetargets:
                    uniquetargets = np.append(uniquetargets, sep)

    for r in rem:
        uniquetargets = uniquetargets[uniquetargets != r]
    uniquetargets = uniquetargets[uniquetargets != "None"]
    uniquetargets = uniquetargets[uniquetargets != ""]
    uniquetargets = uniquetargets.tolist()

    targetmatrix = sparse.lil_matrix((len(targets), len(uniquetargets)), dtype=np.bool)

    for target in targets:
        splittarget = target.split("\t")
        for s in splittarget:
            if target != "None":
                targetmatrix[targets.index(target), uniquetargets.index(s)] = True

    np.savez(DATA + "target_matrix", targetmatrix)
    with open(DATA + "EC_order", "w") as ECfile:
        for target in uniquetargets:
            ECfile.write(target + "\n")

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
