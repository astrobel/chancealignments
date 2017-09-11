import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits as pyfits
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import lomb
import translate as tr
import os, time

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

def translate(value, leftMin, leftMax, rightMin, rightMax):
   # Figure out how 'wide' each range is
   leftSpan = leftMax - leftMin
   rightSpan = rightMax - rightMin

   # Convert the left range into a 0-1 range (float)
   valueScaled = float(value - leftMin) / float(leftSpan)

   # Convert the 0-1 range into a value in the right range.
   return rightMin + (valueScaled * rightSpan)

for i in os.listdir(os.getcwd()):
   if i.endswith("llc.fits"):
      exec("hdulist = pyfits.open('%s')" % i)
      parameters = hdulist[0].header
      kic = parameters['KEPLERID']
      pass
   else:
      continue

exec("lc = np.loadtxt('kic%d_lc.dat')" % kic)
times = lc[:,0]
ampls = lc[:,1]

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
   print ' '
   frequencies, power_spectrum = fft(np.asarray(times), np.asarray(ampls), ofac, hifac)
   hifac = 283 / max(frequencies)
   frequencies, power_spectrum = fft(np.asarray(times), np.asarray(ampls), ofac, hifac)
   power_spectrum = np.concatenate((power_spectrum, [0]))
   power_spectrum = power_spectrum * 4 * np.var(ampls) / ampls.size
   power_spectrum = np.sqrt(power_spectrum)
   power_spectrum *= 1e6 #* ((power_spectrum / (np.median(power_spectrum))) - 1.) # to ppm

   foldfreq = frequencies[power_spectrum.argmax()]

print ' '

tohz = foldfreq * 1e-6
tos = 1/tohz
foldper = tos/86400

binnum = 100
binnum2 = 1000
binsize = foldper / binnum
binsize2 = foldper / binnum2

phasedtimearray = np.zeros(len(times))
finalampls = np.zeros(binnum)
amplcounts = np.zeros(binnum)
finalampls2 = np.zeros(binnum2)
amplcounts2 = np.zeros(binnum2)

for i, val in enumerate(times):
   phasedtime = val % foldper
   newphasedtime = tr.translate(phasedtime, 0, foldper, 0, 1)
   phasedtimearray[i] = newphasedtime
   bindex = (phasedtime - (phasedtime % binsize)) / binsize - 1
   finalampls[bindex] += ampls[i]
   amplcounts[bindex] += 1
   bindex2 = (phasedtime - (phasedtime % binsize2)) / binsize2 - 1
   finalampls2[bindex2] += ampls[i]
   amplcounts2[bindex2] += 1
   
finalampls = np.divide(finalampls, amplcounts)
finalampls2 = np.divide(finalampls2, amplcounts2)

finaltimes = np.histogram(phasedtimearray, bins=binnum-1, range=(0,1))
finaltimes2 = np.histogram(phasedtimearray, bins=binnum2-1, range=(0,1))

plt.figure(1)

plt.plot(finaltimes2[1], finalampls2, 'cs', markersize=5)
plt.plot(finaltimes2[1]+1, finalampls2, 'cs', markersize=5)
plt.plot(finaltimes[1], finalampls, 'ro', markersize=5)
plt.plot(finaltimes[1]+1, finalampls, 'ro', markersize=5)
plt.xlim(0, max(finaltimes2[1]+1))
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
