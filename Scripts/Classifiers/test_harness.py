# A module with functions for testing models
# ------------------------------------------------------------------------------------------------------
# Imports

import numpy as np
from sklearn import metrics
from statistics import mean, pstdev, StatisticsError


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Test the model

def test_model(X_test, model, y_test, ECs):
    pred = predict(X_test, model).round()
    return get_metrics(pred, y_test, ECs)


# ------------------------------------------------------------------------------------------------------
# Produce the metrics of the model
def get_metrics(pred, y_test, ECs):
    top_level = [EC.split(".")[0] for EC in ECs]
    accuracy, precision, sensitivity, specificity = [], [], [], []

    for i in range(1, 8):
        indices = [index for index, value in enumerate(top_level) if value == str(i)]
        act_p = np.where(y_test[:, indices].any(1))[0]
        pred_p = np.where(pred[:, indices].any(1))[0]
        p = np.unique(np.concatenate((act_p, pred_p)))
        accuracy.append(get_accuracy(pred[p, ][:, indices], y_test[p, ][:, indices]))
        print("Accuracy of model group %d.x : %.2f%%" % (i, 100*accuracy[-1]))

        true_positives, true_negatives, false_positives, false_negatives = get_numbers(pred[:, indices], y_test[:, indices])
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
        try:
            sensitivity.append(mean(sensitivities))
            print("Mean sensitivity group %d: (%.2f +/- %.2f)%%" % (i, 100*sensitivity[-1], 100*pstdev(sensitivities)))
            if non_sensitive != 0:
                print("%d EC numbers had no positives in the target matrix, these have been omitted from the sensitivity calculations"
                      % non_sensitive)
        except StatisticsError:
            print("No EC numbers had positives in the target matrix, sensitivity could not be calculated")
        try:
            specificity.append(mean(specificities))
            print("Mean specificity group %d: (%.2f +/- %.2f)%%" % (i, 100*specificity[-1], 100*pstdev(specificities)))
            if non_specific != 0:
                print("%d EC numbers had no negatives in the target matrix, these have been omitted from the sensitivity calculations"
                      " - this is very unlikely" % non_specific)
        except StatisticsError:
            print("No EC numbers had negatives in the target matrix, specificity could not be calculated - this is very unlikely")
        try:
            precision.append(mean(precisions))
            print("Mean precision group %d: (%.2f +/- %.2f)%%" % (i, 100*precision[-1], 100*pstdev(precisions)))
            if non_precise != 0:
                print("%d EC numbers had no positives in the prediction matrix, these have been omitted from the precision calculations"
                      % non_precise)
        except StatisticsError:
            print("No EC numbers had positives in the prediction matrix, precision could not be calculated")
    print("Mean specific accuracy: %.2f%%" % (100*mean(accuracy)))
    print("Mean specific sensitivity: %.2f%%" % (100*mean(sensitivity)))
    print("Mean specific specificity: %.2f%%" % (100*mean(specificity)))
    print("Mean specific precision: %.2f%%" % (100*mean(precision)))
    print("Overall accuracy: %.2f%%" % (100*get_accuracy(pred, y_test)))

    return accuracy, sensitivity, specificity, precision

# -----------------------------------------------------------------------------------------------------
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
    return np.count_nonzero(correct)/correct.size


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
            if y_test[(row, col)] and prediction[(row, col)]:
                t_pos += 1
            elif not y_test[(row, col)] and not prediction[(row, col)]:
                t_neg += 1
            elif y_test[(row, col)] and not prediction[(row, col)]:
                f_neg += 1
            elif not y_test[(row, col)] and prediction[(row, col)]:
                f_pos += 1

        true_positives.append(t_pos)
        true_negatives.append(t_neg)
        false_positives.append(f_pos)
        false_negatives.append(f_neg)

    return true_positives, true_negatives, false_positives, false_negatives


# ------------------------------------------------------------------------------------------------------
# Function to calculate metrics

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

