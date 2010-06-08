#!/usr/bin/env python
from getopt      import getopt
from numpy       import array, mean, std
from scipy.stats import median
import os
import sys

usage = "\nUsage: pywr_validate ix1,[ix2,ix3,..] file1 [file2 file3 ..] ! [file1 [file2 file3 ..]\n\n" +\
        "ix1 .. ixn is a comma delimited list of word indices for which the exact onset matters\n" +\
        "(e.g. YES or NO responses would go here, confidence ratings would not). The script will\n" +\
        "report on onset discrepancies only for these words.\n\n" +\
        "The rest of the arguments are files to compare. The first set of files are those scored\n" +\
        "by the automatic algorithm, these are followed by a separator (!), and then by a list of\n"+\
        "hand scored files. The first automatically scored file is compared to the first manually\n"+\
        "scored file, the second automatically scored file is compared to the second manually\n"+\
        "scored file, and so on."

if len(sys.argv) < 5 or sys.argv.count('!') != 1:
    sys.exit(usage)

separatorIndex  = sys.argv.index('!')
autoScoredFiles = map(os.path.expanduser, sys.argv[2:separatorIndex])
manuScoredFiles = map(os.path.expanduser, sys.argv[separatorIndex+1:])

if len(autoScoredFiles) != len(manuScoredFiles):
    sys.exit(usage)

for file in autoScoredFiles+manuScoredFiles:
    if not os.access(file, os.F_OK):
        sys.exit("The file '%s' doesn't exist." % file)

onsetIndices      = sys.argv[1].split(",")
onsetDistances    = []
numMatched        = 0 
numMatchedUnknown = 0 
numTotal          = 0
numCorrect        = 0

for i in range(len(manuScoredFiles)):
    fd = open(manuScoredFiles[i], 'r')
    manuLines = fd.readlines()
    fd.close()

    fd = open(autoScoredFiles[i], 'r')
    autoLines = fd.readlines()
    fd.close()

    numTotal += max(len(autoLines), len(manuLines))
    numMatchedThisFile = 0

    for manuLine in manuLines:
        # Less autoLines than manuLines (rare).
        if len(autoLines) == 0:
            break

        (manuTime, manuIndex, manuLabel) = manuLine.strip().split("\t")

        # Find closest matching automatic line. 
        bestAutoLine = None
        bestTimeDiff = None
        for (j, autoLine) in enumerate(autoLines):
            (autoTime, autoIndex, autoLabel) = autoLine.strip().split("\t")
            timeDiff = abs(long(manuTime)-long(autoTime))
            if (bestTimeDiff==None or timeDiff < bestTimeDiff) and timeDiff < 250:
                bestAutoLine = j
                bestTimeDiff = timeDiff

        # Matching line?
        if bestTimeDiff == None:
            print "NO MATCH", manuLine, "in", manuScoredFiles[i]
            continue
        numMatchedThisFile += 1

        # Record onset difference.
        (autoTime, autoIndex, autoLabel) = autoLines[bestAutoLine].strip().split("\t")
        if manuIndex in onsetIndices:
            onsetDistances.append(abs(long(manuTime)-long(autoTime)))

        # Is it correct?
        if int(manuIndex) == int(autoIndex):
            numCorrect += 1
        elif int(autoIndex) == -1:
            numMatchedUnknown += 1

        # Delete this auto line from further consideration.
        del autoLines[bestAutoLine]

    numMatched += numMatchedThisFile
    numTotal   += (len(manuLines)-numMatchedThisFile)

    # Unmatched auto lines that are UNKNOWNs are correct.
    for autoLine in autoLines:
        (autoTime, autoIndex, autoLabel) = autoLine.strip().split("\t")
        if int(autoIndex) == -1:
            numCorrect += 1

print "Error rate: (%d/%d) %f" % (numTotal-numCorrect, numTotal, float(numTotal-numCorrect)/float(numTotal))
print "'Unknown' rate: (%d/%d) %f" % (numMatchedUnknown, numMatched, float(numMatchedUnknown)/float(numMatched))
print "Distance mean: %dms std: %dms median: %dms" % (mean(array(onsetDistances)), std(array(onsetDistances)), \
                                                        median(array(onsetDistances)))

