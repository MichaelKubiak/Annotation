#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to identify Swissprot entries with EC numbers
# ------------------------------------------------------------------------------------------------------

from paths import DATA
import re


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    # ------------------------------------------------------------------------------------------------------
    # read in files

    with open(DATA + "enzyme.dat") as efile:
        enzyme = efile.readlines()

    # ------------------------------------------------------------------------------------------------------
    # get a dictionary of EC numbers and the IDs of proteins of that classification

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

    # ------------------------------------------------------------------------------------------------------
    # read in row file from parse_result.py to order list of EC numbers

    with open(DATA + "matrix_rows") as rows:
        order = rows.readlines()

    targetlist = []
    for accession in order:
        accession = accession.strip("\n")
        try:
            targetlist.append(ECs[accession])
        except KeyError:
            targetlist.append("None")

    print(targetlist.count("None")/len(targetlist))

    # ------------------------------------------------------------------------------------------------------
    # output a list of EC numbers which matches the order of the proteins in the matrix

    outstring = ""
    for target in targetlist:
        if target == "None":
            outstring += target + "\n"
        else:
            outstring += "\t".join(target) + "\n"

    with open(DATA + "targets", "w") as targets:
        targets.write(outstring)

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
