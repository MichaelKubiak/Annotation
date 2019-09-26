#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# Script to produce a frequency plot for the number of proteins that are classified as each EC number
# ------------------------------------------------------------------------------------------------------
# Imports

from paths import DATA
import re
import frequency_bars as fb


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to make a dictionary of EC numbers with proteins of that classification

def get_ps_for_EC(enzyme):
    ECs = {}
    current = ""
    for line in enzyme:
        # Set the key for the current dictionary entry
        if line.startswith("ID"):
            current = line.split()[1]
        elif line.startswith("DR"):
            # Add each protein with that classification to the value of the current dictionary entry
            for protein in re.finditer(r"\w+,", line):
                # If the entry already exists, append to a list, otherwise, make the list
                if current in ECs:
                    ECs[current].append(protein.group().split(",")[0])
                else:
                    ECs[current] = [protein.group().split(",")[0]]
    return ECs


# ------------------------------------------------------------------------------------------------------
# Function to produce a logarithmic set of brackets with the number of keys in 'dictionary' that have a number of members in that range

def make_brackets(dictionary, x):
    # Make a dictionary with numbers of items in each value of the input dictionary as the new value
    lengths = {k: len(v) for k, v in dictionary.items()}
    # The first bracket is ECs with 0 proteins
    brackets = [sum(v == 0 for v in lengths.values())]
    # Add 12 further brackets members the previous and current powers of x
    for n in range(1, 13):
        brackets.append(sum(x ** (n-1) <= v < x ** n for v in lengths.values()))
    return brackets


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    # ------------------------------------------------------------------------------------------------------
    # Read in file

    with open(DATA + "enzyme.dat") as efile:
        enzyme = efile.readlines()

    # ------------------------------------------------------------------------------------------------------
    #  Make the frequency plot

    ECs = get_ps_for_EC(enzyme)

    # Set the base for the logarithmic scale
    x = 2

    brackets = make_brackets(ECs, x)

    # Set up labels and positions of the bars
    labels = ["0", "1"] + list("%d - %d" % (x ** n, x ** (n+1)-1) for n in range(1, 12))
    xpos = range(0, 13)

    # Call the function that draws the plot
    fb.freq_bars(brackets, labels, xpos, "linear", "Number of Annotated Proteins", "Number of EC Classifications")


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
