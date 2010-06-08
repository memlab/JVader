#!/usr/bin/env python
import os
import sys
import audio_utils
from optparse import OptionParser

usage = "Usage: pywr_onsets.py [options] file1.wav [file2.wav file3.wav ..]\n\n" +\
        "For each wav file, a fileN.tpa file is generated with the detected onsets."

optParse = OptionParser(usage=usage)
optParse.add_option("--bgFile", help="Path to a .wav file with a recording of typical background noise. " +\
                                     "A default background profile (with little noise) is used if none provided.")
(opts, args) = optParse.parse_args()
if len(args) == 0:
    optParse.error("Please specify at least one .wav file.")

if opts.bgFile != None and not os.access(opts.bgFile, os.F_OK):
    sys.exit("The specified background noise file does not exist.")

files = []
for i in range(len(args)):
    file = os.path.expanduser(args[i])
    if not os.access(file, os.F_OK):
        sys.exit("The file '%s' doesn't exist." % file)
    files.append(file)

for file in files:
    tpaFile = file[:-4] + ".tpa"
    parFile = file[:-4] + ".par"
    if not os.access(tpaFile, os.F_OK) and not os.access(parFile, os.F_OK):
        (dir, filename) = os.path.split(file)
        print "Marking: ", filename[:-4]
        (input, samplingRate) = pywr.loadSpeech(file)
        endpoints = pywr.voicedEndpoints_ns(input, samplingRate, bgFile=opts.bgFile)
   
        fd = open(tpaFile, "w")
        for i in range(len(endpoints)):
            fd.write(str(pywr.ix2ms(endpoints[i][0], samplingRate)) + "\t" +\
                     str(-1) + "\t" +\
                     "[???]\n")
        fd.close()

