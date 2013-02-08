from scipy.io.wavfile import read as wavread
from scipy.io.wavfile import write as wavwrite
from scipy.signal import hann, lfilter, freqz, decimate, butter
from numpy import array, double, amax, absolute, zeros, floor, arange, mean
from numpy import correlate, dot, append, divide, argmax, int16, sqrt, power
from numpy.random import randn
from matplotlib.pyplot import plot, figure, show

def block_process(data, fs, block_len, overlap):
    """
    A generator that slices an input array into overlapping blocks.

    data      the input array as a one dimensional numpy array.
    fs        the sample rate as integer number.
    block_len the length of one block in seconds.
    overlap   the percentage of overlap between blocks.
    """
    block_samples = round(block_len * fs)
    overlap_samples = round(block_samples * overlap)
    shift_samples = block_samples - overlap_samples
    num_blocks = floor((len(data)-overlap_samples) / shift_samples)
    for idx in arange(0, num_blocks):
        samples = data[idx*shift_samples:idx*shift_samples+block_samples]
        yield(array(samples, copy=True), idx*shift_samples)

def levinson_algorithm(r):
    """
    Calculates A and K coefficients of an auto correlation array using
    the levinson-durben method.
    """
    a = zeros(len(r));
    k = zeros(len(r));
    for m in arange(0,len(r)-1):
        alpha = -dot( r[m::-1], append(a[m:0:-1], 1.0) )
        mu = -dot(r[m::-1], append(a[1:m+1], 0.0)) - r[m+1]
        k[m] = -mu / alpha
        a[1:m+2] = append(a[1:m+1], 0.0) + k[m] * append(a[m:0:-1], 1.0)

    a[0] = 1
    return (a,k)

def fundamental_period_estimate(rxx,fs):
    """
    Calculates the fundamental frequency of an auto correlation array.

    rxx   the auto correlation array.
    fs    the sample rate in hertz.
    """
    f_low, f_high = 50 , 250
    f_low_idx = round(fs / f_low)
    f_high_idx = round(fs / f_high)
    period_idx = argmax(rxx[f_high_idx:f_low_idx ]) + f_high_idx
    is_voiced = max(rxx) > 0.20
    return(period_idx, is_voiced)

def vocode(signal, fs, block_len, overlap, order):
    """
    Analyzes a speech signal and synthesizes a vocoded speech signal.

    The speech signal is analyzed using the levinson-durben algorithm
    of the given order. Then, an corresponding output signal is
    synthesized from the levinson-durben coefficients.

    signal     the speech signal as a one dimensional numpy array.
    fs         the sample rate in hertz.
    block_len  the block processing block length in seconds.
    overlap    the block processing block overlap in percent (0..1).
    order      the number of coefficients to use.

    returns a vocoded signal of the same sample rate as the original.
    """

    b, a = butter(1, 200/fs, 'high')
    glottal_lowpass = lambda signal: lfilter(b, a, signal)

    out = zeros(len(signal))
    for block, idx in block_process(signal, fs, block_len, overlap):
        gain_correction = (1-overlap) * 2 # *2 due to hann window
        block *= hann(len(block)) * gain_correction

        rxx = correlate(block, block, mode='full')
        rxx = rxx[len(rxx)/2:]
        period_samples, is_voiced = fundamental_period_estimate(rxx, fs)

        block = preemphasis(block)
        rxx = correlate(block, block, mode='full')
        rxx = rxx[len(rxx)/2:]
        a, k = levinson_algorithm(rxx[:order+1])
        error_power = rms(lfilter(a, (1,), block))

        if is_voiced:
            vocoded = zeros(len(block))
            vocoded[period_samples::period_samples] = 1.0
            vocoded = glottal_lowpass(vocoded)
        else:
            vocoded = randn(len(block))/2

        vocoded = lfilter((error_power,), a, vocoded)
        vocoded *= hann(len(block))
        out[idx:idx+len(block)] += deemphasis(vocoded)
    return out

def preemphasis(signal):
    return lfilter([1, -0.70], 1, signal)

def deemphasis(signal):
    return lfilter([1, 0.70], 1, signal)

def rms(signal):
    return sqrt(mean(power(signal, 2)))

if __name__ == "__main__":
    fs, data = wavread('Mann.wav')
    data = array(data, dtype=double)
    data /= amax(absolute(data))
    data = decimate(data, 4)
    fs = round(fs/4)

    block_len = 0.032
    overlap = 0.5
    order = 16

    out = vocode(data, fs, block_len, overlap, order)

    wavwrite('vocoded.wav', fs, array(out/amax(absolute(out)) * (2**15-1), dtype=int16))

    figure()
    plot(data)
    figure()
    plot(out)
    show()

# ideas:
# use reduce(ola, map(process, array))
# http://stackoverflow.com/questions/6657820/python-convert-an-iterable-to-a-stream
