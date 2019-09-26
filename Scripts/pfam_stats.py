#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# Script to produce a frequency plot for the number of proteins that are classified as hits to each pfam HMM
# ------------------------------------------------------------------------------------------------------
from paths import DATA
from scipy import sparse
import matplotlib.pyplot as plt
import frequency_bars as fb


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to produce a list of numbers of explicit values in each row of a matrix

def get_hits(matrix):
    hits = []
    for row in matrix:
        hits.append(row.nnz)
    print("Max hits on an HMM: %d "%max(hits))
    print("Mean hits on an HMM: %.2f"%(sum(hits)/int(matrix.shape[0])))
    return hits


# ------------------------------------------------------------------------------------------------------
# method to get the indices of explicit values in each column of a matrix and put them in a dictionary entry with their names as the keys

def get_indices(matrix, namelist):
    i = 0
    indices = {}
    for col in matrix:
        indices[namelist[i]] = list(col.indices)
        i += 1
    return indices


# ------------------------------------------------------------------------------------------------------
# method to produce ratios of enzymes to non enzymes for each entry in a dictionary of indices

def get_ratios(indices, targets):
    ratios = []
    nohits = 0
    nonenzyme = []
    enzyme = []
    for hmm, p in indices.items():
        zeros = 0
        nonzeros = 0
        for index in p:
            if targets[index].nnz == 0:
                zeros += 1
            else:
                nonzeros += 1
        if zeros == 0 and nonzeros == 0:
            nohits += 1
        elif zeros == 0:
            enzyme.append(nonzeros)
        #        print("%s has %d hits with EC numbers and no hits without" % (hmm, nonzeros))
        elif nonzeros == 0:
            nonenzyme.append(zeros)
        #        print("%s has %d hits without EC numbers and no hits with" % (hmm, zeros))
        else:
            ratios.append(nonzeros/zeros)
    #        print("%s has %d hits with EC numbers and %d hits with no EC numbers, giving a ratio of %.2f" % (hmm, nonzeros, zeros, (nonzeros/zeros)))
    return nohits, nonenzyme, ratios, enzyme


# ------------------------------------------------------------------------------------------------------
# method to produce a set of brackets for the current plot

def make_brackets(x, enzyme, nonEnzyme, ratios, hits):
    if enzyme is None:
        brackets = [sum(hit == 0 for hit in hits)]
        for n in range(1, 16):
            brackets.append(sum(x ** (n-1) <= hit < x ** n for hit in hits))

        labels = ["0", "1"]+list("%d - %d" % (x ** n, x ** (n+1)-1) for n in range(1, 15))
    else:
        brackets = {"Only Non-enzyme": len(nonEnzyme), "<%d^-3" % x: sum(ratio < x ** -3 for ratio in ratios)}
        for n in range(-2, 6):
            brackets["%d^%d<=r<%d^%d" % (x, n-1, x, n)] = sum(x ** (n-1) <= ratio < x ** n for ratio in ratios)
        brackets["Only Enzymes"] = len(enzyme)
        labels = brackets.keys()
        brackets = brackets.values()
    return brackets, labels


# ------------------------------------------------------------------------------------------------------
# Function to plot a frequency plot

def plot_frequencies(x, yscale, xlabel, ylabel, hits=None, enzyme=None, nonEnzyme=None, ratios=None):

    brackets, labels = make_brackets(x, enzyme, nonEnzyme, ratios, hits)
    xpos = range(0, len(labels))

    fb.freq_bars(brackets, labels, xpos, yscale, xlabel, ylabel)


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    scores_t = sparse.load_npz(DATA + "score_matrix.npz").transpose()

    hits = get_hits(scores_t)
    plot_frequencies(2, "log", "Number of Hits", "Number of Pfam Hidden Markov Models ($log_{10}$)", hits=hits)

    targets = sparse.load_npz(DATA + "target_matrix.npz")

    with open(DATA + "matrix_columns") as pfam_file:
        pfam = pfam_file.readlines()

    pfam = list(map(str.strip, pfam))

    indices = get_indices(scores_t, pfam)

    nohits, nonenzyme, ratios, enzyme = get_ratios(indices, targets)

    plot_frequencies(10, "log", "Ratio of Enzyme to Non-Enzyme Hits", "Number of Pfam Hidden Markov Models ($log_{10}$)",
                     enzyme=enzyme, nonEnzyme=nonenzyme, ratios=ratios)


if __name__ == '__main__':
    main()

