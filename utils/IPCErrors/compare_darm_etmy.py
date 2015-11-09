#!/usr/bin/env python
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

# Modify these to diagnose NDS2 problems
nds_kwargs = {} # {'verbose':True, 'host':'nds.ligo.caltech.edu'

print "Fetching DARM"
darm = TimeSeries.fetch(darm_chan, st, st+dur+1, **nds_kwargs)
print "Fetching ISCINF"
sus = TimeSeries.fetch(sus_chan, st, st+dur+1, **nds_kwargs)
print "Fetching INJ"
inj = TimeSeries.fetch(inj_chan, st, st+dur+1, **nds_kwargs)

# Sample rates all the same
assert(sus.sample_rate.value == darm.sample_rate.value)
assert(inj.sample_rate.value == darm.sample_rate.value)
srate = int(darm.sample_rate.value)

# Hardware inj channel delayed two samples getting to EY SUS
# DARM delayed one sample
# Probably because inj goes from CAL->LSC->SUSEY and DARM just LSC->SUSEY
inj = inj[:-srate]
darm = darm[1:-srate+1]
sus = sus[2:-srate+2]

# SUS channel has opposite sign (ETMY output matrix element probably)
diff = (darm+inj) + sus

# Remove the DARM cal line from the data; it's not recorded in DARM_OUT
import linetools
times, amp, phase = linetools.demodulate(diff, f_cal, 2, 1)
amp = median(amp)
phase = median(phase)
diff_nocal = linetools.subtract_line(diff, amp, f_cal, phase)

glitch_idx = (abs(diff_nocal.value)).argmax()
glitch_t = glitch_idx/float(srate) + st
size = diff_nocal.value[glitch_idx]
print '%.6f amplitude=%.6g' % (glitch_t, size)
print diff_nocal.value[glitch_idx-2:glitch_idx+3]

if False:
    import pylab
    pylab.figure()
    pylab.title('Difference of DARM_OUT and ETMY ISCINF, %u' % (st+pad,))
    pylab.plot(arange(len(diff))/float(srate) - pad, diff_nocal.value)
    pylab.ylabel('Difference (counts)')
    pylab.xlabel('Time (sec)')
    pylab.xlim(-1,2)
    pylab.show()


