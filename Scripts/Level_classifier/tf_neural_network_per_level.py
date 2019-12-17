# A module to build a neural network classifier with tensorflow
# ------------------------------------------------------------------------------------------------------
# Imports

from Regressors import prep
import numpy as np
import tensorflow as tf
from scipy.sparse import vstack
import math


# ------------------------------------------------------------------------------------------------------
# train network - a function to produce the full branched network

def train_network(scores, targets, random_state, n_hidden, n_nodes, activations, output_activation,
                  optimiser, loss, n_epochs, batch_size, ECs):

    X_train, y_train, X_test, y_test = prep.train_test_split_sparse(scores, targets, test_size=0.3, random_state=random_state)
    np.random.seed(random_state)
    tf.random.set_seed(random_state)
    # if only one number of nodes or activation function is given, multiply it by the number of hidden layers so that all hidden layers have those values
    if n_hidden != 1:
        if type(n_nodes) is int:
            n_nodes = [n_nodes] * n_hidden
        if type(activations) is str:
            activations = [activations] * n_hidden
    nn=train_level(ECs, X_train, activations, batch_size, loss, n_epochs, n_hidden, n_nodes, optimiser, output_activation,
                   random_state, y_train)
    # print(nn.evaluate_generator(generate_batch(X_test, y_test, batch_size, n_epochs, random_state, level, ECs), round(len(y_test)/batch_size)))
    return X_test, nn, y_test


# ------------------------------------------------------------------------------------------------------
# train level - a function to produce and train any network of the following - scores are constant, targets are full set - current level produced in batch
#       -   first level  - 8 outputs (None, 1-7)
#       -   second level - x outputs (x = number of second level classifications)
#       -   third level  - ''
#       -   fourth level - ''

def train_level(ECs, X_train, activations, batch_size, loss, n_epochs, n_hidden, n_nodes, optimiser, output_activation,
                random_state, y_train, level=1):
    if level not in range(1, 5):
        raise ValueError(str(level) + " is not a valid value for level please pass a value between 1 and 4")
    nn=tf.keras.Sequential()
    nn.add(tf.keras.layers.Dense(n_nodes[0], activation=activations[0], input_shape=(X_train.shape[1])))
    for i in range(1, n_hidden):
        nn.add(tf.keras.layers.Dense(n_nodes[i], activation=activations[i]))
    nn.add(tf.keras.layers.Dense(y_train.shape[1], output_activation))
    nn.compile(optimizer=optimiser, loss=loss, metrics=["accuracy"])
    nn.fit_generator(generate_batch(X_train, y_train, batch_size, n_epochs, random_state, level, ECs),
                     round(len(y_train)/batch_size), epochs=n_epochs)
    return nn, cat_dict

# ------------------------------------------------------------------------------------------------------
# batch generator - as prior except matrix y

def generate_batch(X_train, y_train, batch_size, n_epochs, random_state, level, ECs):
    # split into each level in turn
    nonone = ECs
    nonone.remove("None")
    prev_level = (0, len(nonone))
    for pos in range(0, level):
        this_level = [EC.split(".")[level] for EC in nonone[prev_level]]
    for i in range(n_epochs):
        for j in range(math.ceil(X_train.shape[0]/batch_size)):
            X_rearr, y_rearr = prep.rearrange(X_train, y_train, random_state+i)
            X_batch = X_rearr[j*batch_size:(j+1)*batch_size]
            y_batch = y_rearr[j*batch_size:(j+1)*batch_size]
