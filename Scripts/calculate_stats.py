#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to produce statistics about the data
# ------------------------------------------------------------------------------------------------------

from paths import DATA
import re

# ------------------------------------------------------------------------------------------------------
# Calculate percentage non enzyme

with open(DATA + "targets") as tfile:
    targets = tfile.readlines()

targets = list(map(str.strip, targets))

print("Non-enzyme: %.2f%%" % (targets.count("None")/len(targets)*100))

# ------------------------------------------------------------------------------------------------------
# Find HMMs not represented in Swissprot

with open(DATA + "hmmresult_full") as resultsfile, open(DATA + "matrix_columns") as pfamfile:
    results, pfam = resultsfile.readlines(), pfamfile.readlines()

unused = []
for line in pfam:
    present = False
    for row in results:
        if not row.startswith("#"):
            splitrow = list(filter(None, re.split(r"\s", row)))
            if splitrow[3] == line:
                present = True

                break
    if not present:
        unused.append(line)
print("Unused HMMs: number = %d %s" % (len(unused), str(unused)))

# ------------------------------------------------------------------------------------------------------
# Find Proteins not represented in Pfam

with open(DATA + "matrix_rows") as spfile:
    swissprot = spfile.readlines()

unmodelled = []
for line in swissprot:
    present = False
    for row in results:
        if not row.startswith("#"):
            splitrow=list(filter(None, re.split(r"\s", row)))
            if splitrow[5] == line:
                present = True
                break
    if not present:
        unmodelled.append(line)
print("Unrepresented proteins: number = %d %s" % (len(unmodelled), str(unmodelled)))
