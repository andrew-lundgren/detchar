from numpy import *
import pylab

# Everything in kg-meter-sec except masses
c = 3.e8
G = 6.67e-11
Msol = 2.e30
Mpc = 3.09e22

# Masses of objects in units of Msol
m1 = 1.4
m2 = 1.4

# Distance to system in Mpc
r = 120.

# CBC signal calculations based on D.Brown's thesis
# http://arxiv.org/abs/0705.1514
# Will denote equations from that
mt = m1+m2
mu = m1*m2/mt # Eq. 2.72

f_isco = 4400./mt # standard rule of thumb

hc_amp = (2./c/c)*mu*((pi*G*Msol*mt)**(2./3.)) # Eq. 2.102
hplus_amp = hc_amp * (G*Msol)/(c*c*Mpc*r)*2. # Eq 2.105
# Note I've just made the system optimally oriented by iota=0
# A factor of f^2/3 in Hz will be added when plotting
# This is the time-domain amplitude at any given frequency
# which is what we want (I think)

# P. Fritschel's tabulation of the maximum drive amplitudes (v1)
# https://dcc.ligo.org/LIGO-T1500484
# I'll just eyeball things from the plots
# Treat them as power laws   y = a x^n

# The ETMY DAC overflow limit
x1, y1 = 70., 1.e-16
x2, y2 = 800., 1.e-19
n_esd = log(y2/y1) / log(x2/x1)
a_esd = y1/(x1**n_esd)
print "ESD coefficients:", a_esd, n_esd

# The Pcal limit
x1, y1 = 100., 2e-16
x2, y2 = 1000., 2e-18
n_pcal = log(y2/y1) / log(x2/x1)
a_pcal = y1/(x1**n_pcal) 
print "PCAL coefficients:", a_pcal, n_pcal

# We'll convert the waveform from strain to meters
# To match Peter's plots
armlength = 4000.

# Output max frequencies 
print "Pcal max f =", exp(log(armlength*hplus_amp/a_pcal)/(n_pcal-2./3.))
print "ESD max f =", exp(log(armlength*hplus_amp/a_esd)/(n_esd-2./3.))

# Do the plotting
freqs = arange(30,2000) 
pylab.figure()
pylab.loglog(freqs, armlength*hplus_amp*(freqs**(2./3.)), c='r', label='Waveform')
pylab.loglog(freqs, a_pcal*(freqs**n_pcal), c='g', label='PCAL limit') 
pylab.loglog(freqs, a_esd*(freqs**n_esd), c='b', label='ETMY ESD limit')
pylab.axvline(f_isco, c='k', ls=':', label='ISCO')
pylab.xlim(30.,2000.)
pylab.ylim(1.e-20,1.e-14)
pylab.xlabel('Frequency (Hz)')
pylab.ylabel('Maximum ETM displacement, peak (m)')
pylab.title('Actuation for a (%.1f, %.1f) CBC injection at %u Mpc' % (m1, m2, r))
pylab.legend(loc='upper right')

pylab.show()

