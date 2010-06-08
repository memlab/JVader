#!/usr/bin/env python

#-- standard lib --#
import array
import cPickle
import math
import optparse
import os
import sys
import wave

#-- external dependencies --#
import numpy
import scipy.signal


def voicedEndpoints_ns(input, samplingRate, smallestSegment=100, useZcCorrection=False, bgFile=None):
    """
    Like voicedEndpoints(..), but first prepends the 
    input with the first 100ms of the .wav file
    specified by 'bgFile'. voicedEndpoints(..) 
    expects the first 100ms to be 'typical' background 
    noise, which it uses to build a template of what 
    background noise in general sounds like. If 'bgFile'
    is None, a default background profile is used.
    """
    if bgFile == None:
        bgFile = "silence.wav"

    (silence, silenceSampleRate) = loadSpeech(bgFile)
    numSilentPoints = math.floor(silenceSampleRate/10.0)

    input = numpy.hstack((silence[0:numSilentPoints], input))
    endpoints = voicedEndpoints(input, samplingRate, smallestSegment, useZcCorrection)

    # Shift endpoints back down by the amount of silence, since
    # it's not really part of the original signal.
    newEndpoints = []
    for (start, end) in endpoints:
        newEndpoints.append((start-numSilentPoints, end-numSilentPoints))

    return newEndpoints

def voicedEndpoints(input, samplingRate, smallestSegment=100, useZcCorrection=False):
    """
    'input' is a 1D array of sound samples at the 
    given sampling rate. Returns a list of 2-tuples 
    containing the beginning and end indices of all
    found utterances. 

    smallestSegment is the length of the smallest
    utterance to keep, in ms. Helps eliminate short
    lived popping / lip smacking.

    If useZcCorrection is True, uses zero-crossing
    statistics in addition to signal energy to refine
    endpoints. This better ensures the voiced segments
    are within the given intervals, but if background
    noise is almost absent and the signal ramps up
    fairly fast in comparison (as is the case in our 
    testing rooms), this has the potential to pad the 
    voiced data with up to 250ms of background on each 
    end; it is thus off by default.
    """
    endpoints = []

    # The algorithm used is similar to the one
    # described in Rabiner and Sambur (1975).
    # We've refined things to detect multiple
    # utterances.
    # We're assuming the SNR ratio is pretty good,
    # which should be true in the testing room.

    # ~10 ms window size.
    windowSize = int(math.floor(samplingRate / 100.0))

    # ~10 ms window step.
    windowStep = int(math.floor(samplingRate / 100.0))
    
    # Smallest segment size in samples.
    smallestSamples = math.ceil((smallestSegment*samplingRate)/1000.0)

    # Smallest segment size in windows.
    smallestWindows = math.ceil(smallestSamples/windowSize)

    # Move the window across the input, count the
    # the number of zero crossings and the energy of 
    # each segment.
    numWindows = int(math.floor((len(input)-windowSize) / windowStep)) + 1
    zeroXings  = numpy.zeros(numWindows)
    energy     = numpy.zeros(numWindows)
    for w in range(numWindows):
        window = input[w*windowStep : w*windowStep + windowSize]
        window = window[window!=0]
        energy[w] = sum(abs(window))
        if useZcCorrection:
            for i in range(len(window)-1):
                if (window[i] < 0 and window[i+1] > 0) or \
                   (window[i] > 0 and window[i+1] < 0):
                    zeroXings[w] += 1

    # Assume first ~100ms is silence.
    numSilentWindows = int(math.floor( ((samplingRate / 10.0)-windowSize) / windowStep )) + 1

    # Used later to refine endpoint estimate, ~250ms in length.
    if useZcCorrection:
        numZcWindows = int(math.floor( ((samplingRate / 4.0)-windowSize) / windowStep )) + 1
    else:
        numZcWindows = 0

    # Zero crossing threshold. If the amount of zero 
    # crossings in a window is more than 2 sd away
    # from the amount in a window of silence, that 
    # window is likely not background noise.
    izct = numpy.mean(zeroXings[0:numSilentWindows]) + 2*numpy.std(zeroXings[0:numSilentWindows])

    # Peak energy.
    imx = numpy.max(energy)

    # Mean energy during silence.
    imn = numpy.mean(energy[0:numSilentWindows])

    # Lower energy threshold. The smaller of 3 percent
    # of the peak energy above silence energy, or four 
    # times the silence energy.
    itl = min(0.03*(imx - imn) + imn, 4*imn)

    # Upper energy thresold.
    itu = 5*itl

    # Look for a window with energy below the lower 
    # threshold.
    p1 = None # Potential starting point.
    n1 = None # Real starting point.
    n2 = None # Real end point.
    for w in range(numWindows):
        # If we don't have a potential starting point and
        # we found a window with energy greater than the
        # lower threshold, mark it as a potential starting point.
        if p1 == None and energy[w] >= itl:
            p1 = w
        # If we have a potential starting point and the current 
        # energy is less than the lower threshold, unmark the
        # the currently marked starting point.
        if p1 != None and energy[w] < itl:
            p1 = None
        # If we haven't found the real starting point, but we have 
        # a potential starting point and we either found a window with 
        # energy greater than the upper threshold or the lower energy
        # threshold has been maintained for more than twice the smallest
        # segment, mark the "potential" starting point as a real starting 
        # point. The criterion of maintaining the lower energy threshold for
        # a period of time without crossing the upper energy threshold is
        # something that did not appear in Rabiner and Sambur (1975),
        # but has been found to be useful with our sound files.
        if n1 == None and p1 != None and (energy[w] >= itu or w - p1 >= 2*smallestWindows):
            n1 = p1
            # As a further refinement, look at the previous 
            # numZcWindows windows. If three or more exceed 
            # the zc threshold, set the starting point to the 
            # window where the zc threshold was first surpassed.
            numSurpassed = 0
            firstWindow  = None
            for wb in range(n1-1, numpy.max(n1-1-numZcWindows,numSilentWindows-1), -1):
                if zeroXings[wb] > izct:
                    numSurpassed += 1
                    firstWindow   = wb
            if numSurpassed >= 3:
                n1 = firstWindow
        # If we have a real starting point, but no ending point,
        # and the current window is below the lower energy threshold,
        # mark the previous window as the ending point.
        if n1 != None and n2 == None and energy[w] < itl:
            n2 = w-1
            # As a further refinement, look at the next
            # numZcWindows windows. If three or more exceed 
            # the zc threshold, set the ending point to the 
            # window where the zc threshold was last surpassed.
            numSurpassed = 0
            lastWindow   = None
            for wb in range(n2+1, min(n2+1+numZcWindows,numWindows)):
                if zeroXings[wb] > izct:
                    numSurpassed += 1
                    lastWindow    = wb
            if numSurpassed >= 3:
                n2 = lastWindow
        # If we have a real starting point, but not ending point,
        # and this is the last frame, mark this as the ending point.
        if n1 != None and n2 == None and w == numWindows-1:
            n2 = w
        # If we found both a starting point and an ending point, add them 
        # to the list and start looking for the next set.
        if n1 != None and n2 != None:
            if n2 - n1 >= smallestWindows:
                pts = (int(math.ceil(n1*windowStep+windowSize/2.0)), int(math.ceil((n2+1)*windowStep-windowSize/2.0)))
                endpoints.append(pts)
            p1 = None
            n1 = None
            n2 = None

    return endpoints
            


def loadSpeech(filename, startPos=0, endPos=None, samplingRate=11025):
    """
    Loads the given wave file, filters out non-speech frequencies,
    downsamples the signal if necessary, and returns a 2-tuple
    with the signal and its sampling rate.
    
    'startPos' is the frame to start reading.
    'endPos' is the frame to stop reading, if None, endPos is taken
    to be the end of the file.
    """
    # Read the file in.
    sizes = "-bh-l"
    f = wave.open(filename, 'rb')
    f.setpos(startPos)
    if endPos == None:
        endPos = f.getnframes()
    origSamplingRate = f.getframerate()
    data = numpy.array(array.array(sizes[f.getsampwidth()], f.readframes(endPos-startPos)) )
    f.close()

    # 2nd-order Butterworth filter.
    low    =  1000
    high   = 16000
    nyq    = origSamplingRate/2.0
    [b, a] = scipy.signal.butter(2, [low/nyq, high/nyq], 'band')
    data   = scipy.signal.lfilter(b, a, data)

    return (decimate(data, origSamplingRate/samplingRate), samplingRate)

def ix2ms(index, samplingRate):
    """
    Given an index into some sound vector and the
    sampling rate for that vector, returns the
    equivalent time rounded to the nearest ms.
    """
    return int(numpy.round((index*1000)/samplingRate))

def ms2ix(time, samplingRate):
    """
    Given time in ms, returns the corresponding
    sample number.
    """
    return int(numpy.round((time*samplingRate)/1000.0))

def decimate(x, q, n=None, ftype='iir', axis=-1):
    """downsample the signal x by an integer factor q, using an order n filter
    
    By default, an order 8 Chebyshev type I filter is used or a 30 point FIR 
    filter with hamming window if ftype is 'fir'.

    (port to python of the GNU Octave function decimate.)

    Inputs:
        x -- the signal to be downsampled (N-dimensional array)
        q -- the downsampling factor
        n -- order of the filter (1 less than the length of the filter for a
             'fir' filter)
        ftype -- type of the filter; can be 'iir' or 'fir'
        axis -- the axis along which the filter should be applied
    
    Outputs:
        y -- the downsampled signal

    This function was borrowed from snippets.dzone.com
    """

    if type(q) != type(1):
        raise Error, "q should be an integer"

    if n is None:
        if ftype == 'fir':
            n = 30
        else:
            n = 8
    if ftype == 'fir':
        b = scipy.signal.firwin(n+1, 1./q, window='hamming')
        y = scipy.signal.lfilter(b, 1., x, axis=axis)
    else:
        (b, a) = scipy.signal.cheby1(n, 0.05, 0.8/q)

        y = scipy.signal.lfilter(b, a, x, axis=axis)

    return y.swapaxes(0,axis)[::q].swapaxes(0,axis)

usage = "Usage: pywr_onsets.py [options] file1.wav [file2.wav file3.wav ..]\n\n" +\
        "For each wav file, a fileN.tpa file is generated with the detected onsets."

optParse = optparse.OptionParser(usage=usage)
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
    (dir, filename) = os.path.split(file)
    print "** marking ", filename + " **"
    (input, samplingRate) = loadSpeech(file)
    endpoints = voicedEndpoints_ns(input, samplingRate, bgFile=opts.bgFile)
    for i in range(len(endpoints)):
        start = str(ix2ms(endpoints[i][0], samplingRate))
        end = str(ix2ms(endpoints[i][1], samplingRate))
        print "(" + start + ", " + end + ")"

