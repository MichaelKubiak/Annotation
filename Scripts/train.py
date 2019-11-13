#!/usr/bin/env python

import Regressors.test_harness as th
import Regressors.train_model as tm
from Regressors import prep
from Regressors.tf_neural_network import train_network

args = tm.arguments("nn", "nn.h5")
proteins, pfam, scores, targets = tm.read_files(args)

ECs, targets = prep.order_ECs(targets)


X_test, network, y_test = train_network(scores, targets, 1, int(args.hidden), int(args.nodes), "tanh", args.epochs, 600, ECs)

network.save(args.path + args.output)
