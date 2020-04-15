# A module to build a neural network classifier with tensorflow
# ------------------------------------------------------------------------------------------------------
# Imports

from Level_classifier import prep
import numpy as np
import tensorflow as tf


# ------------------------------------------------------------------------------------------------------
# train level - train the current level of the network

def train_level(data, targets, level, n_hidden, n_nodes, activations, prev_class="None", random_state=0):

    train, test, unique_targets = prep.train_test_split(targets, prev_class, level, random_state)

    # find number of input features
    n_features = len(train)

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
    nn.add(tf.keras.layers.Dense())
