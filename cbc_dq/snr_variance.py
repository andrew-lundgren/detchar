#! /usr/bin/env python
from gwpy.timeseries import TimeSeries
from numpy import *
import sys

ifo = sys.argv[1]
chan = ifo + ":GDS-CALIB_STRAIN"
st = int(sys.argv[2])
dur = 2048

fmin, fmax = 30., 1024.

chunk, step = 64, 16
fftlen, overlap = 16, 8

data = TimeSeries.fetch(chan, st, st+dur, verbose=True)
srate = data.sample_rate.value

psd_est = data.psd(fftlen, overlap, method='median').value
psd_est[0] = 1.
inv_psd_est = 1. / psd_est
inv_psd_est[0] = 0.

df = 1./float(fftlen)
freqs = df*arange(len(inv_psd_est))
freqs[0] = 1.e-6

integrand = freqs**(-7./3.)*inv_psd_est
integrand[:int(fmin/df)] = 0.
integrand[int(fmax/df):] = 0.

norm = 1./sum(integrand)

result = []
for idx in xrange(0, len(data)-int(srate*chunk), int(srate*step)):
    data1 = data[idx:idx+int(srate*chunk)]
    psd_act = data1.psd(fftlen, overlap).value
    result.append(sqrt(norm*sum(integrand*psd_act*inv_psd_est)))


print result

from matplotlib import use
use('Agg')
import pylab
pylab.figure()
pylab.title("%s, start time %u" % (chan, st))
pylab.axhline(1.)
pylab.plot((step*arange(len(result)) + chunk/2.)/60., result)
pylab.xlabel('Time (min) after start')
pylab.ylabel('Expected SNR standard deviation')
pylab.xlim(0., dur/60.)
pylab.savefig('snr_variance.png')

