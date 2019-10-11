# A module with functions for testing models
# ------------------------------------------------------------------------------------------------------
# Imports

import train_model as tm
import numpy as np
from sklearn import metrics
from statistics import mean, pstdev


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
def test_method(scores, targets, classifier):
    test_scores = []
    for i in range(10):
        print("rand =", i)
        X_test, model, y_test = tm.train_model(scores, targets, i, classifier)
        test_scores.append(test_model(X_test, model, y_test))
    test_scores = np.array(test_scores)
    print("Total mean accuracy:", mean(test_scores[:, 0]))
    print("Total mean sensitivity:", mean(test_scores[:, 1]))
    print("Total mean specificity:", mean(test_scores[:, 2]))
    print("Total mean precision:", mean(test_scores[:, 3]))
    print("Total mean F1 score:", (2*mean(test_scores[:, 1])*mean(test_scores[:, 3])/(mean(test_scores[:, 1]+mean(test_scores[:, 3])))))


# ------------------------------------------------------------------------------------------------------
# Test the model

def test_model(X_test, model, y_test):
    pred = predict(X_test, model)
    accuracy = metrics.accuracy_score(pred, y_test)
    print("Accuracy of model: %.2f%%" % (100*accuracy))
    true_positives, true_negatives, false_positives, false_negatives = get_numbers(pred, y_test)
    sensitivities, specificities, precisions = [], [], []
    non_sensitive, non_specific, non_precise = 0, 0, 0
    for EC in range(len(true_positives)):
        try:
            sensitivities.append(getSensitivity(true_positives[EC], false_negatives[EC]))
        except ZeroDivisionError:
            non_sensitive += 1
        try:
            specificities.append(getSpecificity(true_negatives[EC], false_positives[EC]))
        except ZeroDivisionError:
            non_specific += 1
        try:
            precisions.append(getPrecision(true_positives[EC], false_positives[EC]))
        except ZeroDivisionError:
            non_precise += 1
    sensitivity = mean(sensitivities)
    print("Mean sensitivity: (%.2f +/- %.2f)%%" % (100*sensitivity, 100*pstdev(sensitivities)))
    if non_sensitive != 0:
        print("%d EC numbers had no positives in the target matrix, these have been omitted from the sensitivity calculations"
              % non_sensitive)
    specificity = mean(specificities)
    print("Mean specificity: (%.2f +/- %.2f)%%" % (100*specificity, 100*pstdev(specificities)))
    if non_specific != 0:
        print("%d EC numbers had no negatives in the target matrix, these have been omitted from the sensitivity calculations"
              " - this is very unlikely" % non_specific)
    precision = mean(precisions)
    print("Mean precision: (%.2f +/- %.2f)%%"%(100*precision,100*pstdev(precisions)))
    if non_precise != 0:
        print("%d EC numbers had no positives in the prediction matrix, these have been omitted from the precision calculations"
              % non_precise)
    return accuracy, sensitivity, specificity, precision


# -----------------------------------------------------------------------------------------------------
# Function to predict from the test dataset

def predict(X_test, model):

    # Predict the outcomes of the test data
    prediction = model.predict(X_test)
    return prediction


# # ------------------------------------------------------------------------------------------------------
# # Function to calculate the accuracy of a model - ALREADY IMPLEMENTED IN sklearn.metrics
#
# def get_accuracy(prediction, y_test):
#
#     # Make an equality matrix
#     equality = (prediction == y_test)
#     # Use the equality matrix to produce a matrix showing whether the classifications were correct
#     correct = np.all(equality == True, axis=1)
#     # Determine accuracy from the number of Trues in the matrix
#     return 100*np.count_nonzero(correct)/correct.size


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
    return correct_a/(correct_a + false_b)


# ------------------------------------------------------------------------------------------------------
# Function to calculate sensitivities

def getSensitivity(t_pos, f_neg):
    return calculate_Ratio_True(t_pos, f_neg)


# ------------------------------------------------------------------------------------------------------
# Function to calculate specificities

def getSpecificity(t_neg, f_pos):
    return calculate_Ratio_True(t_neg, f_pos)


# ------------------------------------------------------------------------------------------------------
# Function to calculate precision

def getPrecision(t_pos, f_pos):
    return calculate_Ratio_True(t_pos, f_pos)

