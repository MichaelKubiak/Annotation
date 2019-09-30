# A module with functions to reduce the dataset for testing
# ------------------------------------------------------------------------------------------------------
# imports

from pfam_plots import get_indices, check_targets


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to determine the correct set of HMMs

def choose_hmm(scores, namelist, targets):

    indices = get_indices(scores.transpose(), namelist)

    indices = {k: v for k, v in indices.items() if len(v) == 10}


    newindices = {}
    for k, v in indices.items():
        nonzeros, zeros = check_targets(v, targets)
        if zeros == 0:
            newindices[k] = v
        elif nonzeros/zeros > .8:
            newindices[k] = v

    indexlist = []
    for key in newindices.keys():
        indexlist.append(namelist.index(key))
    test_hmms = scores.transpose()[indexlist]
    return test_hmms.transpose()


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to remove unused EC numbers

def clear_ECs(targets):

    # Create a boolean mask with True where an EC is used

    mask = []
    for EC in targets.transpose():
        if EC.any():
            mask.append(True)
        else:
            mask.append(False)

    # apply the mask to targets to produce the test matrix
    new_targets_t = targets.transpose()[mask]
    return new_targets_t.transpose()



