from gwpy.timeseries import TimeSeries
from gwpy.frequencyseries import FrequencySeries
import numpy as np
from numpy.fft import rfft, irfft
import scipy.signal as sig

def inner_product(data1,data2):
    assert data1.duration==data2.duration
    srate1=data1.sample_rate.value
    srate2=data2.sample_rate.value
    fdata1=rfft(sig.hann(len(data1))*data1.detrend().value)
    fdata2=rfft(sig.hann(len(data2))*data2.detrend().value)
    max_idx=min(len(fdata1),len(fdata2))
    return FrequencySeries(fdata1[:max_idx]*fdata2.conjugate()[:max_idx],
                            df=1./data1.duration)

def avg_freq_bins(fdata,new_df):
    nyquist=(len(fdata)-1)*fdata.df.value
    new_flen=1+int(nyquist/float(new_df))
    binlen=int(new_df/float(fdata.df.value))
    win=sig.hann(2*binlen)
    win/=np.sum(win)
    result=np.zeros(new_flen,dtype=fdata.dtype)
    for idx in range(1,new_flen-1):
        result[idx]=np.sum(win*fdata.value[binlen*(idx-1):binlen*(idx+1)])
    return FrequencySeries(result,df=new_df,channel=fdata.channel,
                            name=fdata.name,epoch=fdata.epoch)

def measure_tf(target, witness, df=1.):
    iprod1=avg_freq_bins(inner_product(target,witness), df)
    iprod2=avg_freq_bins(inner_product(witness,witness), df)
    max_idx=min(len(iprod1),len(iprod2))
    result = iprod1[:max_idx]/iprod2[:max_idx]
    result[0] = 0.
    result[-1] = 0.
    result.name='TF from %s to %s'%(witness.name,target.name)
    return result

def subtract(target, witness, tf):
    assert target.duration == witness.duration
    assert target.epoch == witness.epoch

    srate=int(target.sample_rate.value)
    tlen=len(target)
    flen=1+tlen/2

    tmp=irfft(tf)
    tmp_len=len(tmp)
    tmp=sig.tukey(tmp_len,0.04)*np.roll(tmp, tmp_len/2)
    tmp.resize(tlen) # pad with zeros to length of data
    tf_long=rfft(np.roll(tmp, -tmp_len/2))

    tmp=rfft(witness.detrend().value)
    tmp.resize(flen)

    result=target.detrend().value-irfft(tf_long*tmp)

    pad=ceil(tmp_len/2./srate)
    return TimeSeries(result[int(pad*srate):-int(pad*srate)],
                        sample_rate=target.sample_rate,
                        name='%s minus %s' % (target.name,witness.name),
                        epoch=target.epoch.gps+pad)

