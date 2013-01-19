from scipy.io.wavfile import read as wavread
from scipy.io.wavfile import write as wavwrite
from scipy.signal import hann
from pylab import *

fs, data = wavread('Mann.wav')
block_len = 0.032
overlap = 0.5

def block_process(data, fs, block_len, overlap):
    block_samples = round(block_len*fs)
    overlap_samples = round(block_samples*overlap)
    shift_samples = block_samples - overlap_samples
    num_blocks = floor((len(data)-overlap_samples)/shift_samples)
    for idx in arange(0, num_blocks):
        samples = data[idx*shift_samples:idx*shift_samples+block_samples]
        yield(array(samples, copy=True), idx*shift_samples)

out = zeros(len(data), dtype=int16)
for block, idx in block_process(data, fs, block_len, overlap):
    gain_correction = (1-overlap)*2 # *2 due to hann window
    block *= hann(len(block)) * gain_correction

    # process...
    vocoded = block

    out[idx:idx+len(block)] += vocoded

wavwrite('vocoded.wav', fs, out)

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


if __name__ == "__main__":
    figure()
    plot(data)
    figure()
    plot(out)
    show()

# reduce(ola, map(process, array))
# http://stackoverflow.com/questions/6657820/python-convert-an-iterable-to-a-stream
