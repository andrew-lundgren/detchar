from gwpy.timeseries import TimeSeries
from numpy import *
from numpy.fft import rfft, irfft

# Input: data sample_rate invPSD_length
    # Segment length is just assumed to be length of data

#Steps:

def remove_zeros(arr):
    arr[0] = arr[1]
    arr[-1] = arr[-2]
    for idx in xrange(1,len(arr)-1):
        if arr[idx] == 0.:
            print idx/8.
            arr[idx] = 0.5*(arr[idx-1]+arr[idx+1])
    return arr

# (Over)whiten the data
def whiten(data, sample_rate, whitener_length, overwhiten=True):
    # Estimate PSD
    temp = TimeSeries(data, sample_rate=sample_rate)
    psd = temp.psd(whitener_length, 0.5*whitener_length,method='median').value
    if not overwhiten:
        psd = psd ** 0.5
    # Invert PSD (or ASD, if not overwhitening)
    invpsd = 1./remove_zeros(psd)
    invpsd[0] = 0.
    invpsd[-1] = 0.
    invpsd[:int(16.*whitener_length)] = 0.
    invpsd[int(1025.*whitener_length):] = 0.
    # IFFT the inverse PSD into the time domain
    invpsd_td = irfft(invpsd)
    # Truncate the inverse PSD smoothly, and pad with zeroes    
    trunc_invpsd_td = zeros(len(data),dtype=data.dtype)
    idx = len(invpsd_td)/2
    win = hamming(2*idx)
    trunc_invpsd_td[:idx]=invpsd_td[:idx]*win[idx:]
    trunc_invpsd_td[-idx:]=invpsd_td[-idx:]*win[:idx]
    # FFT the inverse PSD back to frequency domain
    trunc_invpsd = rfft(trunc_invpsd_td)
    # Multiply the data by the inverse FFT
    result = irfft(trunc_invpsd*rfft(data))
    idx = int(sample_rate*whitener_length)
    result[:idx] *= 0.
    result[-idx:] *= 0.
    return psd, invpsd, invpsd_td, trunc_invpsd_td, trunc_invpsd, result

from os.path import expanduser
test_data = TimeSeries.read(expanduser('~/aDQ/data/loud/H1-STRAIN-1118076645-256.hdf'))
psd, invpsd, invpsd_td, trunc_invpsd_td, trunc_invpsd, result = whiten(test_data.value, test_data.sample_rate.value, 8)

import pylab
pylab.plot(arange(len(result))/test_data.sample_rate.value, result)
pylab.show()

