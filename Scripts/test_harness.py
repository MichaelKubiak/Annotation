# A module with functions for testing models
# ------------------------------------------------------------------------------------------------------
# Imports

import numpy as np


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to predict from the test dataset

def predict(X_test, model):

    # Predict the outcomes of the test data
    prediction = model.predict(X_test)
    return prediction


# ------------------------------------------------------------------------------------------------------
# Function to calculate the accuracy of a model

def get_accuracy(prediction, y_test):

    # Make an equality matrix
    equality = (prediction == y_test)
    # Use the equality matrix to produce a matrix showing whether the classifications were correct
    correct = np.all(equality == True, axis=1)
    # Determine accuracy from the number of Trues in the matrix
    return 100*np.count_nonzero(correct)/correct.size


# ------------------------------------------------------------------------------------------------------
# Function to calculate the number of false positives and false negatives

def get_false(prediction, y_test):

    # Make an equality matrix
    equality = (prediction == y_test)
    # For each false in the equality matrix, find whether it was a false positive or false negative
    positive, negative = 0, 0
    for row in range(equality.shape[0]):
        for col in range(equality.shape[1]):
            if not equality[(row, col)]:
                if not prediction[(row, col)]:
                    negative += 1
                elif not y_test[(row, col)]:
                    positive += 1
    return positive, negative



