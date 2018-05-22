from gwpy.timeseries import TimeSeries

ifo='H1'
st_lst=[1160194829,1160194829+16*60]
dur=300

chan_lst=[]
names=[]
for pd in ['AS_A','AS_B','REFL_A','REFL_B']:
    for quad in ['I','Q']:
        for dof in ['PIT','YAW']:
            chan_lst.append('ASC-%s_RF45_%s_%s_OUT_DQ'%(pd,quad,dof))
            names.append(pd+'45'+quad+dof[0])

chan_lst.append('LSC-MOD_RF45_AM_CTRL_OUT_DQ')
names.append('RF45AM')

chan_lst.extend(['LSC-DARM_IN1_DQ','LSC-POP_A_RF45_I_ERR_DQ','LSC-POP_A_RF45_Q_ERR_DQ','LSC-POP_A_RF9_I_ERR_DQ','LSC-POP_A_RF9_Q_ERR_DQ'])
names.extend(['DARM','POP45I','POP45Q','POP9I','POP9Q'])

for st in st_lst:
    for chan,name in zip(chan_lst,names):
        data = TimeSeries.fetch('%s:%s'%(ifo,chan),st,st+dur)
        data.write('%s-%s-%u-%u.hdf'%(ifo,name,st,dur))

