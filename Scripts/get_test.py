#! /usr/bin/env python3
from pathlib import Path
import re

with open(str(Path.home()) + "/Annotation/Data/uniprot_sprot.fasta") as sprotfile:
    sprot = sprotfile.readlines()

i = 0
with open(str(Path.home())+"/Annotation/Data/test_sprot.fasta", "w") as outsprot:
    for line in sprot:
        if line.startswith(">"):
            i += 1
        if i <= 280:
            outsprot.write(line)

print("The swissprot fasta file contains %d sequences" % i)

with open(str(Path.home()) + "/Annotation/Data/Pfam-A.hmm", encoding="utf-8") as pfamfile:
    pfam = pfamfile.readlines()

j = 0
with open(str(Path.home()) + "/Annotation/Data/test_pfam.hmm", "w") as outpfam:
    for line in pfam:
        if line.startswith("//"):
            j += 1
        if j < 90:
            outpfam.write(line)

print("The pfam-A hmm file contains %d hmms" % j)
