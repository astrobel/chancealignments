import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits as pyfits
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import lomb
import smoothing
import translate as tr
import os

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
      r'\usepackage{helvet}',
      r'\usepackage[EULERGREEK]{sansmath}',
      r'\sansmath'
]
mpl.rcParams['axes.formatter.useoffset'] = False

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

print ' '
newphase = raw_input('--> Use highest peak for folding? (y/n) ')
while newphase != "y" and newphase != "n":
   newphase = raw_input('--> Please enter y or n: ')

if newphase == "n":
   while True:
      try:
         foldfreq = float(raw_input('--> Enter folding frequency: (microHertz) '))
         break
      except ValueError:
         print '--> Please enter a float'

elif newphase == "y":
   # second, more intense fourier transform
   print ' '
   print '--> Fourier transform for phase curve <--'
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
   while True:
      try:
         topfreq = float(raw_input('--> Highest frequency to plot in microHertz: '))
         break
      except ValueError:
         print '--> Please enter a float'
   print ' '
   print ' '
   frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
   hifac = topfreq / max(frequencies)
   frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
   power_spectrum = np.concatenate((power_spectrum, [0]))
   power_spectrum = power_spectrum * 4 * np.var(clipped_flux) / clipped_flux.size
   power_spectrum = np.sqrt(power_spectrum)
   power_spectrum *= 1e6 #* ((power_spectrum / (np.median(power_spectrum))) - 1.) # to ppm

   # close-up of interesting part of power spectrum
   plt.figure(1)
   plt.plot(frequencies, power_spectrum, 'r-')
   plt.xlim(0, max(frequencies))
   plt.ylim(ymin = 0)
   plt.xlabel('Frequency ($\mu$Hz)')
   plt.ylabel('Amplitude (ppm)')
   exec("plt.title('%d')" % kic)
   exec("plt.savefig('kic%d_zoomed.png')" % kic)

   foldfreq = frequencies[power_spectrum.argmax()]

# detecting frequency for phase curve and folding
convfactor = (1. / (60 * 60 * 24)) * (10 ** 6)
foldfreq = foldfreq / convfactor
foldper = 1. / foldfreq
clipped_time = clipped_time % foldper

# sort time and flux for binning
sortblend = np.zeros([clipped_time.size, 2])
sortblend = np.array([clipped_time, clipped_flux])
np.reshape(sortblend, (2, clipped_time.size))
sortblend.sort(axis=0)
np.reshape(sortblend, (clipped_time.size, 2))
sortblend = np.transpose(sortblend)
time = sortblend[:,0]
sap_flux = sortblend[:,1]

print ' '
print ' '
print '--> Phase curve <--'
print ' '
smoothtype2 = raw_input('--> Smooth or bin the phase curve? (s/b) ')
while smoothtype2 != "b" and smoothtype2 != "s":
   smoothtype2 = raw_input('--> Please enter s or b: ')

if smoothtype2 == "s": # smoothing

   smoothtype = raw_input('--> Boxcar or Gaussian smoothing? (b/g) ')
   while smoothtype != "b" and smoothtype != "g":
      smoothtype = raw_input('--> Please enter b or g: ')
   while True:
      try:
         kern = int(raw_input('--> Smoothing kernel: '))
         break
      except ValueError:
         print '--> Please enter an integer'
   print ' '

   # smoothing method
   if smoothtype == "b": # boxcar smoothing
      sap_flux2, smth_flux = smoothing.boxsmooth(time, sap_flux, kern)
   elif smoothtype == "g": # gaussian smoothing
      sap_flux2, smth_flux = smoothing.gausssmooth(time, sap_flux, kern)

   # plotting phase curve
   plt.figure(2)
   plt.plot(time, smth_flux, 'ro', markersize=3)
   plt.xlabel('Time mod %f days' % foldper)
   plt.ylabel('Fractional Intensity')
   exec("plt.title('%d')" % kic)
   exec("plt.savefig('kic%d_phase.png')" % kic)

elif smoothtype2 == "b": # binning
   
   while True:
      try:
         binnum = int(raw_input('--> Number of bins: (integer) '))
         break
      except ValueError:
         print '--> Please enter an integer'
   binnum = np.float64(binnum)
   print ' '

   # manual histogram
   binsize = max(time) / binnum
   bindex = 0
   flux_sum = np.zeros(max(time) / binsize)
   flux_num = np.zeros(max(time) / binsize)
   for i, t in enumerate(time):
      bindex = (t - (t % binsize)) / binsize - 1
      flux_sum[bindex] = sap_flux[i] + flux_sum[bindex]
      flux_num[bindex] += 1
   flux_sum = np.divide(flux_sum, flux_num)
   time_binned = np.linspace(0, max(time), binnum)

   time_doubled = np.zeros(time_binned.size)

   for i, val in enumerate(time_binned):
      time_binned[i] = tr.translate(val, 0, foldper, 0, 1)
      time_doubled[i] = time_binned[i] + 1

   binnum2 = np.ceil(binnum / 10.)
   binsize2 = max(time) / binnum2
   bindex2 = 0
   flux_sum2 = np.zeros(max(time) / binsize2)
   flux_num2 = np.zeros(max(time) / binsize2)
   for i, t2 in enumerate(time):
      bindex2 = (t2 - (t2 % binsize2)) / binsize2 - 1
      flux_sum2[bindex2] = sap_flux[i] + flux_sum2[bindex2]
      flux_num2[bindex2] += 1
   time_binned2 = np.linspace(0, max(time), binnum2)
   flux_sum2 = np.divide(flux_sum2, flux_num2)

   time_doubled2 = np.zeros(time_binned2.size)

   for i, val in enumerate(time_binned2):
      time_binned2[i] = tr.translate(val, 0, foldper, 0, 1)
      time_doubled2[i] = time_binned2[i] + 1

   # plotting phase curve
   plt.figure(2)
   plt.plot(time_binned, flux_sum, 'ro', markersize=3)
   plt.plot(time_doubled, flux_sum, 'ro', markersize=3)
   plt.plot(time_binned2, flux_sum2, 'cs')
   plt.plot(time_doubled2, flux_sum2, 'cs')
   plt.xlim(0, max(time_doubled))
   plt.xlabel('Normalised Time mod %f days' % foldper)
   plt.ylabel('Fractional Intensity')
   exec("plt.title('%d')" % kic)
   exec("plt.savefig('kic%d_phase.png')" % kic)

print ' '
show = raw_input('--> Show plot(s) now? (y/n) ')
while show != "y" and show != "n":
   show = raw_input('--> Please enter y or n: ')
print ' '
if show == "y":
   plt.show()
else:
   pass
