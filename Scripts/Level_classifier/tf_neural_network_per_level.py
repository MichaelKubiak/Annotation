# A module to build a neural network classifier with tensorflow
# ------------------------------------------------------------------------------------------------------
# Imports

from Level_classifier import prep
import numpy as np
import tensorflow as tf
import math
from scipy.sparse import vstack


# ------------------------------------------------------------------------------------------------------
# train level - train the current level of the network

def train_level(data, targets, level, n_hidden, n_nodes, activations, n_epochs, batch_size, prev_class="None", random_state=0):

    train, test, unique_targets = prep.train_test_split(targets, 0.3, prev_class, level, random_state)

    # find number of input features
    n_features = data.shape[1]

    # set numpy and tensorflow seeds
    np.random.seed(random_state)
    tf.random.set_seed(random_state)

    # if only one number of nodes or activation function is given, multiply it by the number of hidden layers so that all hidden layers have those values
    if n_hidden != 1:
        if type(n_nodes) is int:
            n_nodes = [n_nodes] * n_hidden
        if type(activations) is str:
            activations = [activations] * n_hidden

    # make a Sequential neural network
    nn = tf.keras.Sequential()
    # add the first hidden layer, which has inputs of the correct shape (the number of input features)
    nn.add(tf.keras.layers.Dense(n_nodes[0], activation=activations[0], input_shape=(n_features,)))
    # create the remaining hidden layers
    for i in range(1, n_hidden):
        nn.add(tf.keras.layers.Dense(n_nodes[i], activation=activations[i]))
    # add the output layer with as many nodes as there are classifications
    nn.add(tf.keras.layers.Dense(len(unique_targets), "sigmoid"))
    # compile the neural network
    nn.compile(optimizer="rmsprop", loss="mse", metrics=["accuracy"])
    # make a generator for the batches to reduce RAM use and fit the network to the generated batches
    nn.fit_generator(generate_batch(data, targets, level, train, batch_size, n_epochs, random_state), round(len(train)/batch_size), epochs=n_epochs)
    # evaluate the network on batches of the test dataset
    print(nn.evaluate_generator(generate_batch(data, targets, test, batch_size, n_epochs, random_state), round(len(test)/batch_size)))


# ------------------------------------------------------------------------------------------------------
# Function to generate a data-target batch of the correct size

def generate_batch(data, targets, level, train, batch_size, epochs=1, random_state=0, keep_duplicate=True):

    # generate a new set of data for each epoch
    for i in range(epochs):
        # rearrange the training dataset and targets
        train_rearr = prep.rearrange(train, random_state+i)
        # generate the correct number of batches for the dataset
        for j in range(math.ceil(len(train)/batch_size)):
            # create the batches
            X_batch = [data[x] for x in train_rearr[j*batch_size:(j+1)*batch_size]]
            # print([targets[x] for x in train_rearr[j*batch_size:(j+1)*batch_size]])
            y_batch = [EC for target in [targets[x] for x in train_rearr[j*batch_size:(j+1)*batch_size]] for EC in [ECs for ECs in target.split("\t")]]
            # for protein in range(len(y_batch)):
            #     if len(y_batch[protein]) > 1:
            #         protein_ECs = y_batch[protein].split("\t")
            #         for EC in range(len(protein_ECs)):
            #             if EC != 0 and keep_duplicate:
            #                 # Duplicate that row of matrix at the bottom of the matrix
            #                 X_batch = vstack(tuple(x for x in X_batch) + tuple(X_batch[0]))
            #                 # add the EC to the end of y_batch
            #                 y_batch.append([protein_ECs[EC]])
            #         # set the original occurrence to the first EC that was present
            #         y_batch[protein] = [protein_ECs[0]]
            # flatten the list into a list of strings, rather than a list of lists
            # y_batch = [y for y_list in y_batch for y in y_list]
            print(y_batch)
            yield X_batch, y_batch



