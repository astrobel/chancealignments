import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits as pyfits
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import lomb
import smoothing
import os

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
      r'\usepackage{helvet}',
      r'\usepackage[EULERGREEK]{sansmath}',
      r'\sansmath'
]

def fft(time, flux, ofac, hifac):
   freq,power, nout, jmax, prob = lomb.fasper(time, flux, ofac, hifac)
   convfactor = (1. / (60 * 60 * 24)) * (10 ** 6)
   uHzfreq = freq * convfactor #11.57, conversion c/d to mHz
   return uHzfreq, power

for i in os.listdir(os.getcwd()):
   if i.endswith("llc.fits"):
      exec("hdulist = pyfits.open('%s')" % i)
      parameters = hdulist[0].header
      kic = parameters['KEPLERID']
      pass
   else:
      continue

# read in light curve
exec("importblend = np.loadtxt('kic%d_lc.dat')" % kic)
clipped_time = importblend[:,0]
clipped_flux = importblend[:,1]

# fourier transform
print ' '
print '--> Fourier transform for power spectrum <--'
print ' '
while True:
   try:
      ofac = int(raw_input('--> Oversampling factor: '))
      break
   except ValueError:
      print '--> Please enter an integer'
while True:
   try:
      hifac = float(raw_input('--> Nyquist range factor: '))
      break
   except ValueError:
      print '--> Please enter a float'
print ' '
frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
hifac = 283 / max(frequencies)
frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
if power_spectrum.size != frequencies.size:
   power_spectrum = np.concatenate((power_spectrum, [0]))
power_spectrum = power_spectrum * 4 * np.var(clipped_flux) / clipped_flux.size
power_spectrum = np.sqrt(power_spectrum)
power_spectrum *= 1e6 #* ((power_spectrum / (np.median(power_spectrum))) - 1.) # to ppm

# logarithmic and linear power spectra
fig1, (pslog, pslin) = plt.subplots(2, 1) 
pslog.plot(frequencies, power_spectrum, 'r-')
pslog.set_xlim(1, max(frequencies))
pslog.set_ylabel('Amplitude (ppm)')
pslog.set_xscale('log')
pslog.set_yscale('log')
exec("pslog.set_title('%d')" % kic)

pslin.plot(frequencies, power_spectrum, 'r-')
pslin.set_xlim(1, max(frequencies))
pslin.set_ylim(ymin = 0)
pslin.set_xlabel('Frequency ($\mu$Hz)')
pslin.set_ylabel('Amplitude (ppm)')

plt.tight_layout()
exec("fig1.savefig('kic%d_spectrum.png')" % kic)

show = raw_input('--> Show plot now? (y/n) ')
while show != "y" and show != "n":
   show = raw_input('--> Please enter y or n: ')
print ' '
if show == "y":
   plt.show()
else:
   pass
