from gwpy.timeseries import TimeSeries
from gwpy.spectrum import Spectrum
from numpy import *

def measure_transfer_function(source, target, tlen, tstride):
    """Measure transfer function from source timeseries to target timeseries. Result is a complex one-sided frequency series. The Nyquist frequency will automatically be set to the smaller of the Nyquist frequencies of source and target.

    source: gwpy TimeSeries
    target: gwpy TimeSeries (same length, can have different sample rate)
    tlen: length of FFT to use in seconds
    tstride: length to step between each FFT (make same as tlen for no overlap)

    Returns: gwpy Spectrum"""

    assert target.epoch == source.epoch
    assert target.duration == source.duration

    srate1 = target.sample_rate.value
    srate2 = source.sample_rate.value

    flen = int(min(srate1*tlen/2+1, srate2*tlen/2+1))

    def fft(data):
        temp = hamming(len(data))*data
        return temp.fft()[:flen]

    crosspower = Spectrum(zeros(flen, dtype=complex128), df=1./tlen)
    sourcepower = Spectrum(zeros(flen, dtype=complex128), df=1./tlen)
    tstarts = range(0,int(target.duration.value-tlen), tstride)
    for tt in tstarts:
        tmp1 = fft(target[int(srate1*tt):int(srate1*(tt+tlen))])
        tmp2 = fft(source[int(srate2*tt):int(srate2*(tt+tlen))])

        crosspower += tmp1*tmp2.conjugate()
        sourcepower += tmp2*tmp2.conjugate()

    transfer = crosspower/sourcepower
    transfer.name = "Transfer function"

    return transfer

