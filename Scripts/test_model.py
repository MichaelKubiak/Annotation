#!/usr/bin/env python

from tensorflow.keras import models
import Regressors.prep as prep
import Regressors.train_model as tm
from Regressors.tf_neural_network import generate_batch
import os
from statistics import mean, pstdev
from paths import DATA
import json

args = tm.arguments("nn", "nn.h5")
proteins, pfam, scores, targets = tm.read_files(args)
ECs, targets = prep.order_ECs(targets)
X_train, X_test, y_train, y_test = prep.train_test_split_sparse(scores, targets, test_size=0.3, random_state=0)

model_file = "/home/mk626/Annotation/Tests/"
evaluation = {}
predictions = []
for directory in os.listdir(model_file):

    nn = models.load_model(model_file + directory)

    # evaluation[directory] = nn.evaluate_generator(generate_batch(X_test, y_test, 600, 1), steps=round(len(y_test)/600))
    # evaluation[directory][0] = float(evaluation[directory][0])
    # evaluation[directory][1] = float(evaluation[directory][1])
    # print(directory, evaluation[directory])

    predictions.append(list(nn.predict_generator(generate_batch(X_test, y_test, 600, 1), steps=round(len(y_test)/600))))

    diff = []
    for i in range(len(predictions[0])):
        print(predictions[0][i])
        print(y_test[i])
        diff.append((predictions[0][i][0]-y_test[i][0])**2)
    print("mean square error: %.2f +/- %.2f" % (mean(diff), pstdev(diff)))

# with open(DATA + "eval_dict", "w") as eval_dict:
#     eval_dict.write(json.dumps(evaluation))
