# ------------------------------------------------------------------------------------------------------
# A module to hold functions for building an SVM regression model that can be used to produce single classifications
# ------------------------------------------------------------------------------------------------------
# Imports

from Regressors import prep
from sklearn.svm import LinearSVR


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Train the SVM

def train_svm(scores, proteins, targets, i):
    # Make learning and test datasets
    X_train, X_test, y_train, y_test, proteins_train, proteins_test = prep.train_test_split_sparse(scores,proteins,targets,test_size=0.3,random_state=i)
    # Make the regressor
    model = LinearSVR(random_state=i)
    model.fit(X_train, y_train)
    return X_test, model, y_test.todense(), proteins_test


# ------------------------------------------------------------------------------------------------------
