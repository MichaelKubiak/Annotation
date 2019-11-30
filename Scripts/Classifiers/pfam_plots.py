#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# Script to produce a frequency plot for the number of proteins that are classified as hits to each pfam HMM
# ------------------------------------------------------------------------------------------------------
# Imports
from paths import DATA
from scipy import sparse
import frequency_bars as fb


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to produce a list of numbers of explicit values in each row of a matrix

def get_hits(matrix):
    hits = []
    # The nnz attribute contains the number of stored values, so the number of hits
    for row in matrix:
        hits.append(row.nnz)
    print("Maximum hits on an HMM: %d " % max(hits))
    print("Mean hits on an HMM: %.2f" % (sum(hits)/int(matrix.shape[0])))
    return hits


# ------------------------------------------------------------------------------------------------------
# Function to get the indices of explicit values in each column of a matrix and put them in a dictionary entry with their names as the keys

def get_indices(matrix, namelist):
    i = 0
    indices = {}
    for col in matrix:
        indices[namelist[i]] = list(col.indices)
        i += 1
    return indices


# ------------------------------------------------------------------------------------------------------
# Function to check whether the rows of a matrix corresponding to a list of values contain values

def check_targets(p, targets):
    zeros = 0
    nonzeros = 0
    for index in p:
        # Check whether the row of targets corresponding to the index has any non zero values, and add one to whichever counter is correct
        if targets[index].nnz == 0:
            zeros += 1
        else:
            nonzeros += 1
    return nonzeros, zeros


# ------------------------------------------------------------------------------------------------------
# Function to produce ratios of enzymes to non enzymes for each entry in a dictionary of indices

def get_ratios(indices, targets):
    ratios = []
    nohits = 0
    nonenzyme = []
    enzyme = []

    # hmm contains a key, p contains the value
    for hmm, p in indices.items():
        nonzeros, zeros = check_targets(p, targets)
        # if both counters are 0, the hmm did not hit any proteins
        if zeros == 0 and nonzeros == 0:
            nohits += 1
        # if zeros is 0, the hmm only hit enzymes
        elif zeros == 0:
            enzyme.append(nonzeros)
        #        print("%s has %d hits with EC numbers and no hits without" % (hmm, nonzeros))
        # if nonzeros is 0, the hmm only hit non-enzymes
        elif nonzeros == 0:
            nonenzyme.append(zeros)
        #        print("%s has %d hits without EC numbers and no hits with" % (hmm, zeros))
        # if both are non-zero, calculate the ratio
        else:
            ratios.append(nonzeros/zeros)
    #        print("%s has %d hits with EC numbers and %d hits with no EC numbers, giving a ratio of %.2f" % (hmm, nonzeros, zeros, (nonzeros/zeros)))
    return nohits, nonenzyme, ratios, enzyme



# ------------------------------------------------------------------------------------------------------
# Function to produce a set of brackets for the current plot

def make_brackets(x, enzyme, nonEnzyme, ratios, hits):
    # if enzyme is None, the plot is for hits per HMM
    if enzyme is None:
        # The first bracket is 0 hits
        brackets = [sum(hit == 0 for hit in hits)]
        # The remaining 15 brackets are logarithmic in x
        for n in range(1, 16):
            brackets.append(sum(x ** (n-1) <= hit < x ** n for hit in hits))
        # create the labels for the bars
        labels = ["0", "1"]+list("%d - %d" % (x ** n, x ** (n+1)-1) for n in range(1, 15))
    else:
        # The brackets are in a dictionary as the easiest way to produce a pair of sets of information
        # The first bracket is non enzyme, followed by less than x ** -3
        brackets = {"Only Non-enzyme": len(nonEnzyme), "<%d^-3" % x: sum(ratio < x ** -3 for ratio in ratios)}
        # Again, the other brackets are logarithmic in x from there up
        for n in range(-2, 6):
            brackets["%d^%d<=r<%d^%d" % (x, n-1, x, n)] = sum(x ** (n-1) <= ratio < x ** n for ratio in ratios)
        # The final bracket is enzymes only
        brackets["Only Enzymes"] = len(enzyme)
        # the dictionary is split into lists
        labels = brackets.keys()
        brackets = brackets.values()
    return brackets, labels


# ------------------------------------------------------------------------------------------------------
# Function to plot a frequency plot

def plot_frequencies(x, yscale, xlabel, ylabel, hits=None, enzyme=None, nonEnzyme=None, ratios=None):
    # Make the brackets
    brackets, labels = make_brackets(x, enzyme, nonEnzyme, ratios, hits)
    # Set positions of the bars
    xpos = range(0, len(labels))
    # Plot the bars
    fb.freq_bars(brackets, labels, xpos, yscale, xlabel, ylabel)


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    # ------------------------------------------------------------------------------------------------------
    # Read in files

    # Load matrix from .npz file
    scores_t = sparse.load_npz(DATA + "score_matrix.npz").transpose()

    with open(DATA + "matrix_columns") as pfam_file:
        pfam = pfam_file.readlines()

    # Strip newlines from the end of list members
    pfam = list(map(str.strip, pfam))

    targets = sparse.load_npz(DATA + "target_matrix.npz")

    # ------------------------------------------------------------------------------------------------------
    # Plot Pfam HMM frequency plot

    hits = get_hits(scores_t)
    plot_frequencies(2, "log", "Number of Hits", "Number of Pfam Hidden \nMarkov Models ($log_{10}$)", hits=hits)

    # ------------------------------------------------------------------------------------------------------
    # Plot Pfam HMM ratio frequency plot

    indices = get_indices(scores_t, pfam)

    nohits, nonenzyme, ratios, enzyme = get_ratios(indices, targets)

    plot_frequencies(10, "log", "Ratio of Enzyme to Non-Enzyme Hits", "Number of Pfam Hidden \nMarkov Models ($log_{10}$)",
                     enzyme=enzyme, nonEnzyme=nonenzyme, ratios=ratios)


# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
