#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to produce statistics about the data
# ------------------------------------------------------------------------------------------------------
# imports

from paths import DATA
import re
import concurrent.futures
import itertools


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Functions to be multithreaded
# ------------------------------------------------------------------------------------------------------
# Find unused hidden markov models

def find_unused(line, results):
    present = False
    line = line.strip("\n")
    for row in results:
        if not row.startswith("#"):
            splitrow = list(filter(None, re.split(r"\s", row)))
            if splitrow[3] == line:
                return
    return line


# ------------------------------------------------------------------------------------------------------
# Find unmodelled proteins

def find_unmodelled(line, results):
    present = False
    line = line.strip("\n")
    for row in results:
        if not row.startswith("#"):
            splitrow = list(filter(None, re.split(r"\s", row)))
            if splitrow[5] == line:
                return
    return line


# ------------------------------------------------------------------------------------------------------
# run batch

def batch(function_name, lines, results):
    output = []
    for line in lines:
        if line is not None:
            output.append(globals()[function_name](line, results))
    return output


# ------------------------------------------------------------------------------------------------------
# itertools recipes grouper

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    # ------------------------------------------------------------------------------------------------------
    # Calculate percentage non enzyme

    with open(DATA + "targets") as tfile:
        targets = tfile.readlines()

    targets = list(map(str.strip, targets))

    print("Non-enzyme: %.2f%%" % (targets.count("None")/len(targets)*100))

    # ------------------------------------------------------------------------------------------------------
    # Find HMMs not represented in Swissprot

    with open(DATA + "hmmresult_56000vs1800") as resultsfile, open(DATA + "matrix_columns") as pfamfile:
        results, pfam = resultsfile.readlines(), pfamfile.readlines()

    unused = []
    executor = concurrent.futures.ProcessPoolExecutor(20)
    futures = [executor.submit(batch, "find_unused", lines, results) for lines in grouper(pfam, 20)]
    concurrent.futures.wait(futures)

    for future in futures:
        unused += future.result()

    while None in unused:
        unused.remove(None)

    print("Unused HMMs: number = %d %s" % (len(unused), str(unused)))

    # ------------------------------------------------------------------------------------------------------
    # Find Proteins not represented by Pfam HMMs

    with open(DATA + "matrix_rows") as spfile:
        swissprot = spfile.readlines()

    unmodelled = []
    futures = [executor.submit(batch, "find_unmodelled", lines, results) for lines in grouper(swissprot, 20)]
    concurrent.futures.wait(futures)
    for future in futures:
        unmodelled += future.result()

    while None in unmodelled:
        unmodelled.remove(None)

    print("Unrepresented proteins: number = %d %s" % (len(unmodelled), str(unmodelled)))

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()
