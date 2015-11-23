from gwpy.timeseries import TimeSeries
from numpy import *
import scipy.signal as sig

cal_lines = {}
cal_lines['H1'] = [35.9, 36.7, 37.3, 331.9] # Is there a PcalX line at 3001.3?
cal_lines['L1'] = [33.7, 34.7, 35.3, 331.3] # Read off a spectrum

# FIXME: Remove this when detrend of TimeSeries is added to gwpy
def detrend(data):
    return TimeSeries(sig.detrend(data), name=data.name,
                        epoch=data.epoch, sample_rate=data.sample_rate)

def window(seglen):
    win = hanning(seglen)
    norm = 2./float(sum(win))
    return norm*win

def make_quadratures(data, f_demod, seglen_secs, stride_secs):
    srate = data.sample_rate.value
    times = arange(len(data), dtype=float64)/float(srate)
    
    data = detrend(data)

    seglen = int(srate*seglen_secs)+1
    stride = int(srate*stride_secs)
    seg_halflen = seglen/2+1

    win = window(seglen)

    cdata = data*cos(2.*pi*f_demod*times)
    sdata = data*sin(2.*pi*f_demod*times)

    samp_times = []
    vals_c = []
    vals_s = []
    for idx in xrange(0, len(data)-seglen, stride):
        samp_times.append(times[idx+seg_halflen])
        vals_c.append(sum(win*cdata[idx:idx+seglen]))
        vals_s.append(sum(win*sdata[idx:idx+seglen]))

    return array(samp_times), array(vals_c), array(vals_s)

def demodulate(data, f_demod, seglen_secs, stride_secs):
    samp_times, vals_c, vals_s = make_quadratures(data, f_demod,
                                                    seglen_secs, stride_secs)

    amp = sqrt(vals_c*vals_c+vals_s*vals_s)
    phase = unwrap(arctan2(-vals_s, vals_c))

    return samp_times, amp, phase

def refine_frequency(data, f_demod, seglen_secs, stride_secs):
    samp_times, vals_c, vals_s = make_quadratures(data, f_demod,
                                                seglen_secs, stride_secs)
    phase = unwrap(arctan2(-vals_s, vals_c))
    a, b = polyfit(samp_times, phase, 1)
    return f_demod + a/(2.*pi)

def subtract_line(data, amp, f_demod, phase):
    srate = data.sample_rate.value
    times = arange(len(data), dtype=float64)/float(srate)
    
    return data - amp*cos(2.*pi*f_demod*times+phase)
 
