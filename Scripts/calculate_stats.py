#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to produce statistics about the data
# ------------------------------------------------------------------------------------------------------

from paths import DATA
import re
import concurrent.futures


# ------------------------------------------------------------------------------------------------------
# Calculate percentage non enzyme

with open(DATA + "targets") as tfile:
    targets = tfile.readlines()

targets = list(map(str.strip, targets))

print("Non-enzyme: %.2f%%" % (targets.count("None")/len(targets)*100))


# ------------------------------------------------------------------------------------------------------
# Function to be multithreaded in next section - finds unused hidden markov models

def find_unused(line, results):
    present = False
    for row in results:
        if not row.startswith("#"):
            splitrow = list(filter(None, re.split(r"\s", row)))
            if splitrow[3] == line:
                present = True

                break
    if not present:
        return line


# ------------------------------------------------------------------------------------------------------
# Find HMMs not represented in Swissprot

with open(DATA + "hmmresult_56000vs1800") as resultsfile, open(DATA + "matrix_columns") as pfamfile:
    results, pfam = resultsfile.readlines(), pfamfile.readlines()

unused = []
executor = concurrent.futures.ProcessPoolExecutor(20)
futures = [executor.submit(find_unused, line, results) for line in pfam]
concurrent.futures.wait(futures)

for future in futures:
    unused.append(future.result())

print("Unused HMMs: number = %d %s" % (len(unused), str(unused)))


# ------------------------------------------------------------------------------------------------------
# Function to be multithreaded in next section - finds unmodelled proteins

def find_unmodelled(line, results):
    present = False
    for row in results:
        if not row.startswith("#"):
            splitrow=list(filter(None, re.split(r"\s", row)))
            if splitrow[5] == line:
                present = True
                break
    if not present:
        return line


# ------------------------------------------------------------------------------------------------------
# Find Proteins not represented by Pfam HMMs

with open(DATA + "matrix_rows") as spfile:
    swissprot = spfile.readlines()

unmodelled = []
futures = [executor.submit(find_unmodelled, line, results) for line in swissprot]
concurrent.futures.wait(futures)
for future in futures:
    unmodelled.appen(future.result())
print("Unrepresented proteins: number = %d %s" % (len(unmodelled), str(unmodelled)))
