# A module with functions to reduce the dataset for testing
# ------------------------------------------------------------------------------------------------------
# imports

from Regressors.pfam_plots import get_indices, check_targets
from Regressors.prep import remove_non_family


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to determine the correct set of HMMs

def choose_hmm(scores, pfam, targets):

    # get indices of explicit values in each column of the score matrix
    indices = get_indices(scores.transpose(), pfam)

    # use only hmms that have exactly 10 (arbitrary number that reduces the matrix size significantly) proteins
    indices = {k: v for k, v in indices.items() if len(v) == 10}

    # Use only hmms with a high proportion of enzyme hits (for more effective testing)
    newindices = {}
    for k, v in indices.items():
        nonzeros, zeros = check_targets(v, targets)
        if zeros == 0:
            newindices[k] = v
        elif nonzeros/zeros > .8:
            newindices[k] = v

    # Make a list of the indices of pfam hmms that will be used
    indexlist = []
    for key in newindices.keys():
        indexlist.append(pfam.index(key))
    # Use that list to reduce the score matrix and list of pfam names down to the test set
    test_hmms = scores.transpose()[indexlist]
    pfam = [pfam[index] for index in indexlist]
    return test_hmms.transpose(), pfam


# ------------------------------------------------------------------------------------------------------
# Create test dataset

def create_test(scores, proteins, pfam, targets):
    scores, pfam = choose_hmm(scores, pfam, targets)
    test_dataset, proteins, pfam, test_targets = remove_non_family(scores, proteins, pfam, targets)

    return test_dataset, proteins, pfam, test_targets

