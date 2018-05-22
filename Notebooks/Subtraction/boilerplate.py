from gwpy.timeseries import TimeSeries
from gwpy.plotter import FrequencySeriesPlot
from numpy import *
from numpy.fft import rfft, irfft
import scipy.signal as sig
import matplotlib.pyplot as plt

ifo='H1'
st1=1160194829
st2=1160194829+16*60
dur=300

all_names=['DARM','POP9Q']
for pd in ['AS_A','AS_B','REFL_A','REFL_B']:
    for quad in ['I','Q']:
        for dof in ['PIT','YAW']:
            all_names.append(pd+'45'+quad+dof[0])

names=['DARM','POP9Q','AS_B45IP','AS_A45QP','AS_B45QY','REFL_B45IY']

data_lst=[]
for name in names:
    print 'Get', name
    data_lst.append(TimeSeries.read('%s-%s-%u-%u.hdf'%(ifo,name,st2,dur)).detrend())

