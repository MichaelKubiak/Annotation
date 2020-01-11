# A module to build a neural network classifier with tensorflow
# ------------------------------------------------------------------------------------------------------
# Imports

from Regressors import prep
import numpy as np
import tensorflow as tf
from scipy.sparse import vstack
import math


# ------------------------------------------------------------------------------------------------------
# train level - train the current level of the network

def train_level(prev_class, level, hidden, nodes, activation):
