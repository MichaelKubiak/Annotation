#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to build a neural network classifier with tensorflow
# ------------------------------------------------------------------------------------------------------
# Imports

import prep
from test_data import create_test
import joblib
import test_harness as th
import numpy as np
from statistics import mean
import train_model as tm
import tensorflow as tf


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Train a neural network classifier

def train_network(scores, targets, i, n_hidden, n_nodes, activations):

    # Make learning and test datasets
    X_train, X_test, y_train, y_test = prep.train_test_split_sparse(scores, targets, test_size=0.3, random_state=i)
    n_features = X_train.shape[1]
    n_classes = y_train.shape[1]
    np.random.seed(i)
    tf.set_random_seed(i)
    if n_hidden != 1:
        if type(n_nodes) is int:
            n_nodes = [n_nodes] * n_hidden
        if type(activations) is str:
            activations = [activations] * n_hidden
    layers = [n_features]
    nn = tf.keras.Sequential()
    for i in range(n_hidden):
        nn.add(tf.keras.layers.Dense(n_nodes[i], activation=activations[i], input_shape=(layers[-1],)))
        layers.append(n_nodes[i])
    nn.add(tf.keras.layers.Dense(n_classes, activation="sigmoid"))
    nn.compile(optimizer="nadam", loss="binary_crossentropy")
    nn.fit(X_train, y_train, epochs=1000, verbose=1)
    print(th.test_model(X_test, nn, y_test))
    print("done")
    # g = tf.Graph()
    # with g.as_default():
    #     tf.set_random_seed(i)
    #     tf_x = tf.placeholder(dtype=tf.float32, shape=(None, n_features))
    #     tf_y = tf.placeholder(dtype=tf.bool, shape=(None, n_classes))
    #
    #     if n_hidden != 1:
    #         if len(n_nodes) == 1:
    #             n_nodes = [n_nodes] * n_hidden
    #         if len(activations) == 1:
    #             activations = [activations] * n_hidden
    #
    #     layers = [tf_x]
    #     for i in range(n_hidden):
    #         layers.append(tf.layers.dense(inputs=layers[-1], units=n_nodes[i], activation=activations[i]))
    #     layers.append(tf.layers.dense(inputs=layers[-1], units=n_classes, activation=None))
    #
    #     predictions = tf.sparse_softmax(layers[-1], axis=1)
    #
    #     with g.as_default():
    #         cost = tf.compat.v2.losses.BinaryCrossentropy()


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    args = tm.arguments("Neural Network", "neural_network.clf")

    # ------------------------------------------------------------------------------------------------------
    # Read files

    proteins, pfam, ECs, scores, targets = tm.read_files(args)

    # ------------------------------------------------------------------------------------------------------
    # make a test dataset

    scores, targets = create_test(scores, pfam, targets)

    # ------------------------------------------------------------------------------------------------------
    # Remove proteins with no pfam hits - nothing happens with test set

    print("Percentage empty rows in target matrix before pruning: %.2f%%" % (prep.get_empty(targets)))

    scores, targets = prep.remove_non_family(scores, targets)

    print("Percentage empty rows in target matrix after pruning: %.2f%%" % (prep.get_empty(targets)))

    # ------------------------------------------------------------------------------------------------------
    # Remove non-enzyme proteins down to a limit

    # limit = 0.2

    # scores, targets = prep.remove_non_enzyme(scores, targets, limit)

    # print("Percentage empty rows in target matrix after removal of empty rows down to %.2f: %.2f%%" % (limit, prep.get_empty(targets)))
    # ------------------------------------------------------------------------------------------------------
    # test method
    test_scores = []
    for i in range(10):
        print("rand =", i)
        X_test, forest, y_test = train_network(scores, targets, i)
        test_scores.append(th.test_model(X_test, forest, y_test))
    test_scores = np.array(test_scores)
    print("Total mean accuracy:", mean(test_scores[:, 0]))
    print("Total mean sensitivity:", mean(test_scores[:, 2]))
    print("Total mean specificity:", mean(test_scores[:, 3]))
    print("Total mean precision:", mean(test_scores[:, 1]))
    print("Total mean F1 score:", (2*mean(test_scores[:, 2])*mean(test_scores[:, 1])/(mean(test_scores[:, 2] + mean(test_scores[:, 1])))))
    # Output the classifier as a pickle using joblib
    X_test, network, y_test = train_network(scores, targets, 1)
    joblib.dump(network, args.path + args.output)

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
