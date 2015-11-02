from gwpy.timeseries import TimeSeries
from gwpy.table.lsctables import SnglBurstTable
from glue import datafind
from numpy import *
from scipy.interpolate import UnivariateSpline
import sys

ifo = sys.argv[1]
chan = sys.argv[2]
st = int(sys.argv[3])
dur = int(sys.argv[4])
et = st + dur

primary_channel = '%s:%s' % (ifo, chan)
vco_channel = ifo + ":SYS-TIMING_C_FO_A_PORT_11_SLAVE_CFC_FREQUENCY_5"

# VCO is sampled in the middle of the second, reported one second later
conn = datafind.GWDataFindHTTPConnection()
cache = conn.find_frame_urls(ifo[0], '%s_R' % ifo, st, et, urltype='file')
vco = TimeSeries.read(cache, vco_channel, st, et+2)
vco = vco[8::16].value
vco_major_MHz = int(median(vco)/1.e6)
vco -= 1e6*vco_major_MHz
vco_interp = UnivariateSpline(arange(len(vco))-0.5, vco)

trigs = SnglBurstTable.fetch(primary_channel, 'omicron', st, et,
                                filt=lambda x: st <= x.get_peak() < et)

result = []
for trg in trigs:
    f_vco = vco_interp(float(trg.get_peak())-st)
    result.append((f_vco, trg.peak_frequency, trg.snr))
save('%s-%s_versus_vco-%u-%u.npy' % (ifo, chan, st, dur), result)

