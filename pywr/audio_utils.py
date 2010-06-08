from endpoints    import voicedEndpoints_ns
from array        import array as pyarray
from numpy        import argmax, array, hstack, inf, max, round
from scipy.signal import butter, cheby1, firwin, lfilter
import cPickle
import os
import wave

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

