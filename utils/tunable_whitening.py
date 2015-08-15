#!/usr/bin/env python
from gwpy.timeseries import TimeSeries
import scipy.signal as sig
from numpy import *
from numpy.fft import irfft

def detrend(data):
    return TimeSeries(sig.detrend(data), epoch=data.epoch,
                            sample_rate=data.sample_rate)

# Measure ASD, given some data
def measure_asd(data, fftlensec, overlap_frac=None, method='median-mean'):
    # Detrend
    dtdata = detrend(data)
    
    if method == 'median-mean':
        print 'Warning: Using default overlap for median-mean method'
        overlap = None
    elif overlap_frac:
        overlap = overlap_frac*fftlensec
    else:
        overlap = None

    asd = dtdata.asd(fftlength=fftlensec, overlap=overlap, method=method)

    return asd

def build_whitener(data, fftlensec, overlap_frac='None', method='median-mean'):
    asd = measure_asd(data, fftlensec, overlap_frac, method)
    invasd = 1./asd
    invasd[0] = 0.
    invasd[-1] = 0. # Might as well set Nyquist to zero
    
    if True: # Get rid of the high-frequency yuckiness
        idx = int(3000./asd.df.value)
        invasd[idx+1:] = invasd[idx]

    return invasd

def apply_whitening(data, whitener):
    fdata = data.fft()
    wdata = TimeSeries(irfft(fdata * whitener), epoch=data.epoch,
                            sample_rate=data.sample_rate)
    print type(wdata)
    return wdata

def impulse_like(data):
    temp = zeros(len(data))
    temp[len(temp)/2] = 1.
    return TimeSeries(temp, epoch=data.epoch, sample_rate=data.sample_rate)

def step_like(data):
    temp = zeros(len(data))
    temp[len(temp)/2:] = 1.
    return TimeSeries(temp, epoch=data.epoch, sample_rate=data.sample_rate)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 6:
        print "Args: chan t_PSD dur_PSD seglen t_glitch"
        exit()
    chan = sys.argv[1]
    t1_psd = int(sys.argv[2])
    dur_psd = int(sys.argv[3])
    seglen = int(sys.argv[4])
    tt = float(sys.argv[5])

    st = int(tt) - seglen/2
    
    data_for_psd = TimeSeries.fetch(chan, t1_psd, t1_psd+dur_psd)
    data = TimeSeries.fetch(chan, st, st+seglen)
    invasd = build_whitener(data_for_psd, seglen,  method='median-mean')
    data = step_like(data)
    temp = apply_whitening(data, invasd)
    #plot = data.plot()
    plot = temp.plot()
    plot.show()

