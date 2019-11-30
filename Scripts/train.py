#!/usr/bin/env python

import Regressors.test_harness as th
import Regressors.train_model as tm
from Regressors import prep
from Regressors.tf_neural_network import train_network
from Regressors.test_data import create_test

args = tm.arguments("nn", "nn.h5")
proteins, pfam, scores, targets = tm.read_files(args)

scores, proteins, pfam, targets = create_test(scores, proteins, pfam, targets)
ECs, targets = prep.order_ECs(targets)


X_test, network, y_test = train_network(scores, targets, 1, int(args.hidden), int(args.nodes), "tanh", args.epochs, 20, ECs)



network.save(args.path + args.output)
