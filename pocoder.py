from scipy.io.wavfile import read as wavread
from scipy.io.wavfile import write as wavwrite
from scipy.signal import hann, lfilter, freqz, decimate
from pylab import *

def block_process(data, fs, block_len, overlap):
    block_samples = round(block_len*fs)
    overlap_samples = round(block_samples*overlap)
    shift_samples = block_samples - overlap_samples
    num_blocks = floor((len(data)-overlap_samples)/shift_samples)
    for idx in arange(0, num_blocks):
        samples = data[idx*shift_samples:idx*shift_samples+block_samples]
        yield(array(samples, copy=True), idx*shift_samples)

def levinson_algorithm(r):
    a = zeros(len(r));
    k = zeros(len(r));
    for m in arange(0,len(r)-1):
        alpha = -dot( flipud(r[0:m+1]), append(flipud(a[1:m+1]),1.0) );
        mu = -dot(flipud(r[0:m+1]),append(a[1:m+1],0.0)) - r[m+1];
        k[m] = -divide(mu,alpha);
        a[1:m+2] = append(a[1:m+1],0.0) + k[m]*append(flipud(a[1:m+1]),1.0);

    a[0] = 1;
    return (a,k)

def fundamental_period_estimate(rxx,fs):
    f_low,f_high = 80 , 250
    f_low_idx = round(fs/f_low)
    f_high_idx = round(fs/f_high)
    period_idx = argmax(rxx[f_high_idx:f_low_idx ])+f_high_idx
    is_voiced = max(rxx)>0.25
    return(period_idx,is_voiced)

def vocode(signal, fs, block_len, overlap, order):
    out = zeros(len(signal))
    for block, idx in block_process(signal, fs, block_len, overlap):
        gain_correction = (1-overlap)*2 # *2 due to hann window
        block *= hann(len(block)) * gain_correction
        rxx = correlate(block,block,mode='full')
        rxx = rxx[len(rxx)/2:]

        a, k = levinson_algorithm(rxx[:order+1])
        period_samples,is_voiced = fundamental_period_estimate(rxx,fs)

        if is_voiced:
            vocoded = zeros(len(block))
            vocoded[period_samples::period_samples] =  1.0
        else:
            vocoded = randn(len(block))

        vocoded = lfilter(rxx[0:1], a, vocoded)
        vocoded *= hann(len(block))
        out[idx:idx+len(block)] += vocoded
    return out

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

    wavwrite('vocoded.wav', fs, array(out/amax(absolute(out)) * 2**14, dtype=int16))

    figure()
    plot(data)
    figure()
    plot(out)
    show()

# reduce(ola, map(process, array))
# http://stackoverflow.com/questions/6657820/python-convert-an-iterable-to-a-stream
