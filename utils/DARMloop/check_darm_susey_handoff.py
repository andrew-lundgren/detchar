from gwpy.timeseries import TimeSeries
from numpy import *
import sys

ifo = 'L1'
f_cal = {'L1': 33.7, 'H1': 37.3}[ifo]
pad = 8
st = int(sys.argv[1]) - pad
dur = 2*pad

darm_chan = ifo + ":LSC-DARM_OUT_DQ"
sus_chan = ifo + ":SUS-ETMY_L3_ISCINF_L_IN1_DQ"
inj_chan = ifo + ":CAL-INJ_HARDWARE_OUT_DQ"

nds_kwargs = {'verbose':False, 'host':'nds.ligo-la.caltech.edu'}

print "Fetching DARM"
darm = TimeSeries.fetch(darm_chan, st, st+dur+1, **nds_kwargs)
print "Fetching ISCINF"
sus = TimeSeries.fetch(sus_chan, st, st+dur+1, **nds_kwargs)
print "Fetching INJ"
inj = TimeSeries.fetch(inj_chan, st, st+dur+1, **nds_kwargs)
assert(sus.sample_rate.value == darm.sample_rate.value)
assert(inj.sample_rate.value == darm.sample_rate.value)
srate = int(darm.sample_rate.value)


inj = inj[:-srate]
darm = darm[1:-srate+1]
sus = sus[2:-srate+2]

# SUS channel has opposite sign (ETMY output matrix element probably)
diff = (darm+inj) + sus

# Remove the cal line from the SUS data (it's not recorded in DARM_OUT)
import linetools
times, amp, phase = linetools.demodulate(diff, f_cal, 2, 1)
amp = median(amp)
phase = median(phase)
print amp, phase
diff_nocal = linetools.subtract_line(diff, amp, f_cal, phase)


glitch_t = (abs(diff_nocal.value)).argmax()/float(srate) + st
size = max(abs(diff_nocal.value))
print '%.6f %.6g' % (glitch_t, size)

import pylab
pylab.figure()
pylab.title('Difference of DARM_OUT and ETMY ISCINF, %u' % (st+pad,))
pylab.plot(arange(len(diff))/float(srate) - pad, diff_nocal.value)
pylab.ylabel('Difference (counts)')
pylab.xlabel('Time (sec)')
pylab.xlim(-1,2)
pylab.show()

if False:
    from IPython import embed
    embed()

