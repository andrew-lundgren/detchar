from gwpy.timeseries import TimeSeries
from numpy import *
import sys

ifo = 'L1'
f_cal = 33.7
pad = 8
st = int(sys.argv[1]) - pad
dur = 2*pad

nds_kwargs = {'verbose':False}

omc_chan = ifo + ":OMC-DCPD_SUM_OUT_DQ"
darm_chan = ifo + ":LSC-DARM_IN1_DQ"

omc = TimeSeries.fetch(omc_chan, st, st+dur+1, **nds_kwargs)
darm = TimeSeries.fetch(darm_chan, st, st+dur+1, **nds_kwargs)

xarr = omc.value - mean(omc.value)
yarr = darm.value - mean(darm.value)

a, b = polyfit(xarr, yarr, 1)

diff = yarr - (a*xarr+b)

print max(abs(diff))

if True:
    import pylab
    pylab.plot(diff)
    pylab.show()

