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

def get_numbers(prediction, y_test):

    # Make an equality matrix
    equality = (prediction == y_test)

    # For each EC number, determine the number of true positives, true negatives, false positives and false negatives
    true_positives, true_negatives, false_positives, false_negatives = [], [], [], []
    for col in range(equality.shape[1]):
        t_pos, t_neg, f_pos, f_neg = 0, 0, 0, 0
        for row in range(equality.shape[0]):
            if equality[(row, col)] and prediction[(row, col)]:
                t_pos += 1
            elif equality[(row, col)] and not prediction[(row, col)]:
                t_neg += 1
            elif prediction[(row, col)]:
                f_neg += 1
            else:
                f_pos += 1

        true_positives.append(t_pos)
        true_negatives.append(t_neg)
        false_positives.append(f_pos)
        false_negatives.append(f_neg)

    return true_positives, true_negatives, false_positives, false_negatives


# ------------------------------------------------------------------------------------------------------
# Function to calculate sensitivities or specificities

def calculate_Ratio_True(correct_a, false_b):
    try:
        return correct_a/(correct_a + false_b)
    except ZeroDivisionError:
        return 1


# ------------------------------------------------------------------------------------------------------
# Function to calculate sensitivities

def getSensitivity(t_pos, f_neg):
    return calculate_Ratio_True(t_pos, f_neg)


# ------------------------------------------------------------------------------------------------------
# Function to calculate specificities

def getSpecificity(t_neg, f_pos):
    return calculate_Ratio_True(t_neg, f_pos)





