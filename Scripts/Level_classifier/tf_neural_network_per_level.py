# A module to build a neural network classifier with tensorflow
# ------------------------------------------------------------------------------------------------------
# Imports

from Level_classifier import prep
import numpy as np
import tensorflow as tf


# ------------------------------------------------------------------------------------------------------
# train level - train the current level of the network

def train_level(data, targets, level, n_hidden, n_nodes, activation, prev_class="None", random_state=0):

    X_train, X_test, y_train, y_test = prep.train_test_split(data, targets, prev_class, level, random_state)
