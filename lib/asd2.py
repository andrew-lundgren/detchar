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
    tmp=sig.hann(tmp_len)*np.roll(tmp, tmp_len/2)
    tmp.resize(len(witness)) # pad with zeros to length of witness
    tf_long=rfft(np.roll(tmp, -tmp_len/2))

    tmp=rfft(witness.astype(np.float64).value)*tf_long
    tmp.resize(flen)

    result=target.value-irfft((srate/float(witness.sample_rate.value))*tmp)

    pad=int(np.ceil(0.5/tf.df.value))
    print pad, srate
    return TimeSeries(result[int(pad*srate):-int(pad*srate)],
                        sample_rate=target.sample_rate,
                        name='%s minus %s' % (target.name,witness.name),
                        epoch=target.epoch.gps+pad)

def bp_chunk(ts, new_srate, f_low, f_high):
    srate=int(ts.sample_rate.value)
    bp=sig.firwin(4*srate,[f_low, f_high],nyq=srate/2.,window='hann',pass_zero=False)
    bp.resize(len(ts))
    tmp=(2.*new_srate/float(srate))*abs(rfft(bp))*rfft(ts)
    padidx=4*int(new_srate) # Sample rate after irfft is twice desired
    return TimeSeries(irfft(tmp[:1+int(ts.duration.value*new_srate)])[padidx:-padidx:2],
                        sample_rate=new_srate,epoch=ts.epoch.gps+2)

def _split_thirds(somearr):
    skip=len(somearr)/3
    return np.concatenate([somearr[:skip],somearr[-skip:]])

def _smooth_mid(myarr):
    skip=len(myarr)/3
    xarr=np.linspace(-1.,1.,3*skip)
    a,b,c=np.polyfit(_split_thirds(xarr),
                     _split_thirds(myarr),
                     2)
    myarr[skip:-skip]=a*xarr[skip:-skip]*xarr[skip:-skip]+b*xarr[skip:-skip]+c
    
def clean_tf(mytf,f1,f2):
    assert f2 > f1
    df=mytf.df.value
    idx1,idx2 = int(f1/df), int(f2/df)
    width=idx2-idx1
    idx0,idx3 = idx1-width, idx2+width
    
    result=mytf.copy()
    
    _smooth_mid(result.real[idx0:idx3])
    _smooth_mid(result.imag[idx0:idx3])
    
    return result

