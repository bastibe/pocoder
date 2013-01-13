from scipy.io.wavfile import read as wavread
from pylab import *

fs, data = wavread('Mann.wav')
block_len = 0.032
overlap = 0.75

def block_process(data, fs, block_len, overlap):
    block_samples = round(block_len*fs)
    overlap_samples = round(block_samples*overlap)
    shift_samples = block_samples - overlap_samples
    num_blocks = floor((len(data)-overlap_samples)/shift_samples)
    for idx in arange(0, num_blocks):
        samples = data[idx*shift_samples:idx*shift_samples+block_len]
        yield(samples, idx)

for block, idx in block_process(data, fs, block_len, overlap):
    print(idx)
    consumer.eat(block)

reduce(ola, map(process, array))

# http://stackoverflow.com/questions/6657820/python-convert-an-iterable-to-a-stream
