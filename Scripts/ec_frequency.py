#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# Script to produce a frequency plot for the number of proteins that are classified as each EC number
# ------------------------------------------------------------------------------------------------------
from paths import DATA
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# function to make a dictionary of EC numbers with proteins of that classification

def get_ps_for_EC(enzyme):
    ECs = {}
    current = ""
    for line in enzyme:
        if line.startswith("ID"):
            current = line.split()[1]
        elif line.startswith("DR"):
            for protein in re.finditer(r"\w+,", line):
                if current in ECs:
                    ECs[current].append(protein.group().split(",")[0])
                else:
                    ECs[current] = [protein.group().split(",")[0]]
    return ECs


# ------------------------------------------------------------------------------------------------------
# function to produce a logarithmic set of brackets with the number of keys in 'dictionary' that have a number of members in that range

def get_brackets(dictionary, x):

    lengths = {k: len(v) for k, v in dictionary.items()}
    brackets = [sum(v == 0 for v in lengths.values())]
    for n in range(1, 13):
        brackets.append(sum(x ** (n-1) <= v < x ** n for v in lengths.values()))
    return brackets


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    with open(DATA + "enzyme.dat") as efile:
        enzyme = efile.readlines()

    ECs = get_ps_for_EC(enzyme)

    x = 2

    brackets = get_brackets(ECs, x)

    t = "bar"

    if t == "line":
        a = [0] + list(x ** n for n in range(12))


        plt.plot(a, brackets)
        plt.xscale('symlog', basex=x)

        plt.gca().xaxis.set_major_formatter(ScalarFormatter())
        plt.show()
    elif t == "bar":
        labels = ["0", "1"] + list("%d - %d" % (x ** n, x ** (n+1)-1) for n in range(1, 12))
        xpos = range(0, 13)

        plt.bar(xpos, brackets, edgecolor="black")
        plt.xticks(xpos, labels, rotation=10)
        plt.xlabel("Number of Annotated Proteins")
        plt.ylabel("Number of EC Classifications")
        plt.show()

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
