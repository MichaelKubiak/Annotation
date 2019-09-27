# Annotation

## Modules

#### paths.py

A variable containing the path to the data folder

#### frequency_bars.py

A function to draw consistent bar charts

#### get_test.py

Functions to produce a test dataset

#### prep.py

Functions to:

- reduce the data by removing proteins that do not hit any HMMs

- remove a portion of the non enzyme proteins, as they will be less important to the learning process


## Scripts

### bash

#### get_data

Downloads the data files:

- Pfam-A.hmm          - containing the hidden markov models from Pfam A
                         
- enzyme.data         - containing the EC numbers and which proteins are members of those groups
                         
- uniprot_sprot.fasta - containing the sequences of the proteins in swissprot

### python

#### parse_result.py

Uses the output of hmmsearch (profile against sequence database) to generate a sparse matrix of hit scores between proteins and pfam hmms

#### identify.py

Produces a sparse, boolean matrix with Trues where a protein is annotated as having a particular EC number

#### calculate_stats.py

Calculates:

1. how much of the swissprot database does not have EC numbers
2. how many (and which) Pfam HMMs are not seen in swissprot
3. how many (and which) swissprot proteins do not have a family in pfam

#### ec_frequency.py

Plots a frequency bar chart for number of proteins per EC number, with logarithmic bar widths

#### pfam_plots.py

Plots frequency bar charts for:

- Number of hits per HMM

- Ratio of enzyme to non-enzyme hits per HMM

