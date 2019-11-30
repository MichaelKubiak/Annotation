# A module to build a neural network classifier with tensorflow
# ------------------------------------------------------------------------------------------------------
# Imports

from Regressors import prep
from Regressors import test_harness as th
import numpy as np
import tensorflow as tf
from scipy.sparse import vstack
import math


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to train a neural network classifier

def train_network(scores, targets, random_state, n_hidden, n_nodes, activations, n_epochs, batch_size, ECs):

    # Make learning and test datasets
    X_train, X_test, y_train, y_test = prep.train_test_split_sparse(scores, targets, test_size=0.3, random_state=random_state)
    # find number of input features
    n_features = X_train.shape[1]
    # n_classes = max(list(map(max, y_train)))
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
    # add the output layer with a single node (regression) with rectified linear unit activation function because the output is bounded
    nn.add(tf.keras.layers.Dense(1, "relu"))
    # compile the neural network with the chosen optimiser and loss function
    nn.compile(optimizer="rmsprop", loss="mse", metrics=["accuracy"])
    # make a generator for the batches to reduce RAM use and fit the network to the generated batches
    nn.fit_generator(generate_batch(X_train, y_train, batch_size, n_epochs, random_state), round(len(y_train)/batch_size), epochs=n_epochs)
    # evaluate the network on batches of the test dataset
    print(nn.evaluate_generator(generate_batch(X_test, y_test, batch_size, n_epochs, random_state), round(len(y_test)/batch_size)))
    #print(th.test_model(X_test, nn, y_test, ECs))
    return X_test, nn, y_test


# ------------------------------------------------------------------------------------------------------
# Function to generate a data-target batch of the correct size

def generate_batch(X_train, y_train, batch_size, epochs=1, random_state=0, keep_duplicate=True):
    # generate a new set of data for each epoch
    for i in range(epochs):
        # generate the correct number of batches for the dataset
        for j in range(math.ceil(X_train.shape[0]/batch_size)):
            # rearrange the training dataset and targets
            X_rearr, y_rearr = prep.rearrange(X_train, y_train, random_state+i)
            # create the batches
            X_batch = X_rearr[j*batch_size:(j+1)*batch_size]
            y_batch = y_rearr[j*batch_size:(j+1)*batch_size]
            # split each protein with multiple ECs into separate lines
            for protein in range(len(y_batch)):
                if len(y_batch[protein]) > 1:
                    for EC in range(len(y_batch[protein])):
                        if EC != 0 and keep_duplicate:
                            # Duplicate that row of matrix at the bottom of the matrix
                            X_batch = vstack((X_batch, X_batch[protein]))
                            # add the EC to the end of y_batch
                            y_batch.append([y_batch[protein][EC]])
                    # set the original occurence to the first EC that was present
                    y_batch[protein] = [y_batch[protein][0]]
            # flatten the list into a list of strings, rather than a list of lists
            y_batch = [y for y_list in y_batch for y in y_list]
            yield X_batch, y_batch

