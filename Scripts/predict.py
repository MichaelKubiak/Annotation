#! /usr/bin/env python3
# ------------------------------------------------------------------------------------------------------
# A script to predict from a previously pickled classifier
# ------------------------------------------------------------------------------------------------------
# Imports

import argparse
from paths import DATA
import joblib


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(description="A script to predict EC number of a protein vector based on a previously pickled classifier")
    parser.add_argument("-c", "--classifier", help="Name of classifier file")
    args = parser.parse_args()
    classifier = joblib.load(DATA + args.classifier)
    prediction = classifier.predict()


# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

# Don't run if imported
if __name__ == '__main__':
    main()
