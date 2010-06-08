from endpoints    import voicedEndpoints_ns
from mfcc         import mfcc
from paths        import DATA_INSTALL_PATH
from array        import array as pyarray
from numpy        import argmax, array, hstack, inf, max, round
from scipy.signal import butter, cheby1, firwin, lfilter
import cPickle
import os
import wave

PYWR_GARBAGE_LABEL = 'PYWR_UNKNOWN'

def soundSamplingRate(filename):
    """
    Returns the sampling rate of the given file (only
    .wav is currently supported).
    """
    f = wave.open(filename, 'rb')
    samplingRate = f.getframerate()
    f.close()
    return samplingRate

def soundNumFrames(filename):
    """
    Returns the current number of frames in the given file
    (only .wav is currently supported). This is mostly 
    useful to mark segments of a file being written to in 
    real time.
    """
    f = wave.open(filename, 'rb')
    numFrames = f.getnframes()
    f.close()
    return numFrames    

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
    data = array( pyarray(sizes[f.getsampwidth()], f.readframes(endPos-startPos)) )
    f.close()

    # 2nd-order Butterworth filter.
    low    =  1000
    high   = 16000
    nyq    = origSamplingRate/2.0
    [b, a] = butter(2, [low/nyq, high/nyq], 'band')
    data   = lfilter(b, a, data)

    return (decimate(data, origSamplingRate/samplingRate), samplingRate)

def train(input, samplingRate, label, verbose=False, bgFile=None):
    """
    Create, fit, and return a new word model using 
    the given input sequence for training.
    A separate model should be trained for each
    word to recognize.

    'input' is a sound vector containing utterances 
    (separated by short periods of silence) of the 
    word to train.

    'label' specifies the classification label for 
    the given class of utterances.

    If verbose is set to True, details of the 
    training will be printed.

    'bgFile' is the path of a .wav file with a
    recording of typical background noise. See
    voicedEndpoints_ns(..).
    """
    endpoints = voicedEndpoints_ns(input, samplingRate, bgFile=bgFile)
    if len(endpoints) == 0:
        print "Warning: training data for", label, "has no vocalizations!"
        return None
    elif verbose:
        print "Found %d utterances in %s" % (len(endpoints), label)

    observations = []
    for (start, end) in endpoints:
        observations.append(mfcc(input[start:end], samplingRate))

    numStates   = 3
    numMixtures = 1
    hmm = GausMixHMM()
    hmm.setLabel(label)
    hmm.setVerbose(verbose)
    hmm.initParams(numStates, numMixtures, observations)
    hmm.baumWelch(observations)

    return hmm

def trainDirs(dirs, words, verbose=False, bgFile=None):
    """
    Looks in all 'dirs' for training data for
    each word (expected in dir/word.wav), 
    concatenates all of the available data, 
    and fits a model to it. A list of 
    word model objects is returned.

    'dirs' is a list of directory paths.
    
    'words' is a list of words, one model is returned
    per word.

    If verbose is set to True, details of the 
    training will be printed.

    'bgFile' is the path of a .wav file with a
    recording of typical background noise. See
    voicedEndpoints_ns(..).

    NOTE: For simplicity, we assume all .wav files are 
    recorded with the same sampling rate. 
    """
    models = []
    for word in words:
        # Grab all available data for this word.
        data = array([])
        for dir in dirs:
            wavFile = os.path.join(dir, word) + ".wav"
            if os.access(wavFile, os.F_OK):
                (dirData, samplingRate) = loadSpeech(wavFile)
                data = hstack((data, dirData))

        # Fit parameters to this word.
        if verbose:
            print "Training: ", word
        model = train(data, samplingRate, word, verbose, bgFile)
        if model != None:
            models.append(model)

    return models

def saveModels(models, dir):
    """
    Saves the given models in 'dir'.

    'models' is a list of word model objects.

    'dir' is the directory to save them in. There's one
    file generated for each model.
    """
    for model in models:
        modelFile = os.path.join(dir, model.getLabel()) + ".hmm"
        fd = open(modelFile, "w")
        cPickle.dump(model, fd) 
        fd.close()

def loadModels(path):
    """
    Returns a list of model objects, stored at the 
    specified path, that can later be given to 
    classify(..).

    'path' is the directory which contains the *.hmm 
    files to load. All *.hmm files found in this 
    directory are loaded. 
    """
    models = []
    for filename in os.listdir(path):
        if filename[-4:].lower() == ".hmm":
            fd    = file(os.path.join(path, filename), "r")
            model = cPickle.load(fd)
            fd.close()
            models.append(model)
   
    return models

def loadGarbageModel():
    """
    Returns the default garbage model object.
    """
    fd = file(os.path.join(DATA_INSTALL_PATH, PYWR_GARBAGE_LABEL+".hmm"), "r")
    garbageModel = cPickle.load(fd)
    fd.close()

    return garbageModel

def classify(input, samplingRate, models, garbageModel, verbose=False, bgFile=None, realtime=False):
    """
    'input' is a sound vector containing utterances
    (separated by short periods of silence) that need
    to be classified. 'models' is a list of word
    models as returned by loadModels(..) and should
    contain a model for each word that may appear
    in the input sequence.

    'bgFile' is the path of a .wav file with a
    recording of typical background noise. See
    voicedEndpoints_ns(..).

    If 'realtime' is set to True, at most one
    label is returned. We loop through the 
    detected utterances backwards and return
    the label of the first valid (non-garbage)
    utterance.

    Returns a list of endpoints and classification 
    labels for all found utterances.
    """
    labels     = []
    numUnknown = 0
    endpoints  = voicedEndpoints_ns(input, samplingRate, bgFile=bgFile)
    if verbose:
        print "Found %d utterances" % len(endpoints)

    for (start, end) in reversed(endpoints):
        obs          = mfcc(input[start:end], samplingRate)
        allStateLike = []

        # Compute the likelihood of the observation under each valid
        # word model.
        for model in models:
            (foo, stateLike) = model.viterbi(obs)
            allStateLike.append(stateLike)

        # Compute likelihood under the garbage model.
        (foo, garbageStateLike) = garbageModel.viterbi(obs)

        # If the best matched valid word model is a better match than
        # the garbage model, go with the best valid word.
        if max(allStateLike) > garbageStateLike:
            label = models[argmax(allStateLike)].getLabel()
            # If realtime, return right away.
            if realtime:
                return ([(start, end)], [label])
        # Otherwise label as garbage.
        else:
            label = garbageModel.getLabel()
            numUnknown += 1

        labels.append(label)
    labels.reverse()   
 
    if verbose:
        print "Marked as unknown: ", numUnknown

    # If we are in real-time mode and we got here, that means
    # there are no valid labels. Don't return anything (as
    # as opposed to e.g. a list of 'UNKNOWN' labels).
    if realtime:
        return ([], [])

    return (endpoints, labels)

def ix2ms(index, samplingRate):
    """
    Given an index into some sound vector and the
    sampling rate for that vector, returns the
    equivalent time rounded to the nearest ms.
    """
    return int(round((index*1000)/samplingRate))

def ms2ix(time, samplingRate):
    """
    Given time in ms, returns the corresponding
    sample number.
    """
    return int(round((time*samplingRate)/1000.0))

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
        b = firwin(n+1, 1./q, window='hamming')
        y = lfilter(b, 1., x, axis=axis)
    else:
        (b, a) = cheby1(n, 0.05, 0.8/q)

        y = lfilter(b, a, x, axis=axis)

    return y.swapaxes(0,axis)[::q].swapaxes(0,axis)

