from numpy        import arange, array, dot, empty, exp, fix, float64, log10, ones, pi, real, sqrt, zeros
from numpy.fft    import fft
from scipy.signal import hamming, triang, lfilter

def mfcc(input, samplingRate):
    """
    Same as mfccE, but with energy term removed.
    """
    return mfccE(input, samplingRate)[:,1:]

def mfccE(input, samplingRate):
    """
    Computes the Mel Frequency Cepstral Coefficients for 
    the given input. 'input' is a 1D array of samples
    at the given sampling rate.

    Returns a matrix where each row is a window of the 
    original signal and each column is a mfcc for that 
    window.
    """
    # The general algorithm for computing mfcc is as follows:
    # 
    # 1. Move a window across the original signal. For each window,
    # 2.    take the DFT
    # 3.    convert the result to the mel scale: create an overlapping 
    #       filter bank of triangular windows (several linearly spaced
    #       ones followed by several log spaced ones, see below) and
    #       convert each resulting magnitude to the log scale
    # 4.    take the discrete cosine transform of the result to pick 
    #       out the independent features
    # 
    # The implementation and some assumptions are inspired by 
    # the Auditory Toolbox for MATLAB by M. Slaney. 

    # Parameters for the filter bank.
    lowestFrequency = 133.3333
    linearFilters   = 13
    linearSpacing   = 66.666666
    logFilters      = 27
    logSpacing      = 1.0711703
    totalFilters    = linearFilters + logFilters

    # Parameters for each window.
    fftSize         = 512
    windowSize      = 256
    windowStep      = 100
    cepstralCoeffs  = 13

    # Lower, center, and upper frequencies for each band of the
    # filter bank.
    freqs                                  = zeros((1,totalFilters+2))
    freqs[0][0:linearFilters]              = lowestFrequency + arange(linearFilters)*linearSpacing
    freqs[0][linearFilters:totalFilters+2] = freqs[0][linearFilters-1] * logSpacing**(arange(logFilters+2)+1)

    lower  = freqs[0][0:totalFilters]
    center = freqs[0][1:totalFilters+1]
    upper  = freqs[0][2:totalFilters+2]

    # Pre-compute triangular weights for each band.
    fftFreqs     = (samplingRate*arange(fftSize))/fftSize
    trianHeight  = 2.0/(upper-lower)
    trianWeights = zeros((totalFilters, fftSize))
    for fi in range(totalFilters):
        trianWeights[fi][:] = (array(fftFreqs > lower[fi],int) & array(fftFreqs <= center[fi],int))* \
                                trianHeight[fi] * (fftFreqs-lower[fi])/(center[fi]-lower[fi])    + \
                              (array(fftFreqs > center[fi],int) & array(fftFreqs < upper[fi],int)) * \
                                trianHeight[fi] * (upper[fi]-fftFreqs)/(upper[fi]-center[fi])

    # Pre-emphasis filter.
    preEmphasized = lfilter(array([1, -0.97]), array([1]), input)

    numWindows = int(fix((len(input)-windowSize)/windowStep))
    ceps       = zeros((numWindows, cepstralCoeffs))

    # For each window.
    for wi in range(numWindows):
        # Indices of the window in the original signal.
        first   = wi*windowStep
        last    = first+windowSize

        # Zero-pad and multiply by a hamming window.
        fftData = zeros((1,fftSize))
        fftData[0][0:windowSize] = preEmphasized[first:last] * hamming(last-first)
    
        # step 2 from above
        fftMag  = abs(fft(fftData)).transpose()

        # steps 3 and 4 from above
        ceps[wi][:] = dct(log10(dot(trianWeights, fftMag)).transpose())[0][0:cepstralCoeffs]

    return ceps

def dct(x,axis=-1):
    """
    Discrete cosine transform based on the FFT.

    For even-length signals it uses an N-point FFT
    For odd-length signals it uses a 2N-point FFT.
    """
    # This implementation is taken from the SciPy "sandbox".
    # Code in the sandbox is considered experimental and is not
    # installed by default (although this implementation seems
    # solid), so we put it here for convenience. We also make
    # a small modification (see below).

    n = len(x.shape)
    N = x.shape[axis]
    even = (N%2 == 0)
    slices = [None]*4
    for k in range(4):
        slices[k] = []
        for j in range(n):
            slices[k].append(slice(None))
    if even:
        xtilde = 0.0*x
        slices[0][axis] = slice(None,N/2)
        slices[1][axis] = slice(None,None,2)
        slices[2][axis] = slice(N/2,None)
        slices[3][axis] = slice(N,None,-2)
    else:
        newshape = list(x.shape)
        newshape[axis] = 2*N
        xtilde = empty(newshape,float64)
        slices[0][axis] = slice(None,N)
        slices[2][axis] = slice(N,None)
        slices[3][axis] = slice(None,None,-1)

    for k in range(4):
        slices[k] = tuple(slices[k])
    xtilde[slices[0]] = x[slices[1]]
    xtilde[slices[2]] = x[slices[3]]
    Xt = fft(xtilde,axis=axis)
    
    pk = exp(-1j*pi*arange(N)/(2*N))

    # The following two lines do not appear in the
    # SciPy dct implementation. We do this to make 
    # the final matrix orthogonal.
    pk[0] /= sqrt(2.)
    pk *= sqrt(2./N)

    newshape = ones(n)
    newshape[axis] = N
    pk.shape = newshape

    if not even:
        pk /= 2;
        Xt = Xt[slices[0]]
    
    return real(Xt*pk)

