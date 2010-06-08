#!/usr/bin/env python
import pywr
import os
import sys
import cPickle
from numpy    import array, hstack
from optparse import OptionParser

usage = "Usage: pywr_train.py [options] training_directory1 [training_directory2 ..]\n\n" +\
        "Looks in training_directory1 for a set of wav files to train on. There should\n" +\
        "be one word per wav file repeated several times with a short silence between\n" +\
        "each repetition. training_directory2 and onward are searched for each of the\n" +\
        "files found in training_directory1, and each word model is trained using all\n" +\
        "all of the available data. A .hmm file is generated in training_directory1\n" +\
        "for each file."

optParse = OptionParser(usage=usage)
optParse.add_option("--bgFile", help="Path to a .wav file with a recording of typical background noise. " +\
                                     "A default background profile (with little noise) is used if none provided.")
(opts, args) = optParse.parse_args()
if len(args) == 0:
    optParse.error("Please specify at least one training directory.")

# Background noise file should exist, if given.
if opts.bgFile != None and not os.access(opts.bgFile, os.F_OK):
    sys.exit("The specified background noise file does not exist.")

# At least training_directory1 should exist.
trainingDirs = map(lambda x: os.path.expanduser(x), args)
if not os.access(trainingDirs[0], os.F_OK):
    sys.exit(trainingDirs[0] + " does not exist.")

# Look in training_directory1 for files to train on.
words = []
for filename in os.listdir(trainingDirs[0]):
    if filename[-4:].lower() == ".wav":
        modelFile = os.path.join(trainingDirs[0], filename[:-4]) + ".hmm"
        # If we haven't yet fit a model to this word, 
        # add the word to the list.
        if not os.access(modelFile, os.F_OK):
            words.append(filename[:-4])

# Train on the words found above.
models = pywr.trainDirs(trainingDirs, words, True, opts.bgFile)

# Save models in training_directory1.
pywr.saveModels(models, trainingDirs[0])

