import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits as pyfits
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import smoothing
import os

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
      r'\usepackage{helvet}',
      r'\usepackage[EULERGREEK]{sansmath}',
      r'\sansmath'
]

print ' '
print '--> Smoothing and clipping light curve <--'
print ' '
smoothtype = raw_input('--> Boxcar or Gaussian smoothing? (b/g) ')
while smoothtype != "b" and smoothtype != "g":
   smoothtype = raw_input('--> Please enter b or g: ')
while True:
   try:
      kern = int(raw_input('--> Smoothing kernel: '))
      break
   except ValueError:
      print '--> Please enter an integer'

kic = 0
quarters = []

for i in os.listdir(os.getcwd()):
   if i.endswith("llc.fits"):
     
      exec("hdulist = pyfits.open('%s')" % i)

      parameters = hdulist[0].header
      
      if kic == 0:
         kic = parameters['KEPLERID']
      else:
         pass

      q = parameters['QUARTER']
      quarters.append(q)

      table = hdulist[1].data

      sap_flux = table['SAP_FLUX']
      time = table['TIME']

      hdulist.close()

      # creating a blend array to remove NaNs
      blend = np.array([time, sap_flux])
      blend = np.transpose(blend)
      blend2 = np.ma.compress_rows(np.ma.fix_invalid(blend))

      exec("time_%d = blend2[:,0]" % q)
      exec("sap_flux_%d = blend2[:,1]" % q)

      if smoothtype == "b": # boxcar smoothing
         exec("sap_flux2_%d, smth_flux_%d = smoothing.boxsmooth(time_%d, sap_flux_%d, kern)" % (q, q, q, q))
      elif smoothtype == "g": # gaussian smoothing
         exec("sap_flux2_%d, smth_flux_%d = smoothing.gausssmooth(time_%d, sap_flux_%d, kern)" % (q, q, q, q))

      continue
   else:
      continue

quarters.sort()

### PLOTTING ###

# unsmoothed data
plt.figure(1)
for i, q in enumerate(quarters):
   if q % 5 == 0:
      exec("plt.plot(time_%d, sap_flux_%d, 'ro', time_%d, smth_flux_%d, 'c-', markersize=3)" % (q, q, q, q))
   if q % 5 == 1:
      exec("plt.plot(time_%d, sap_flux_%d, 'yo', time_%d, smth_flux_%d, 'c-', markersize=3)" % (q, q, q, q))
   if q % 5 == 2:
      exec("plt.plot(time_%d, sap_flux_%d, 'go', time_%d, smth_flux_%d, 'c-', markersize=3)" % (q, q, q, q))
   if q % 5 == 3:
      exec("plt.plot(time_%d, sap_flux_%d, 'bo', time_%d, smth_flux_%d, 'c-', markersize=3)" % (q, q, q, q))
   if q % 5 == 4:
      exec("plt.plot(time_%d, sap_flux_%d, 'mo', time_%d, smth_flux_%d, 'c-', markersize=3)" % (q, q, q, q))
plt.xlabel('Time (d)')
plt.ylabel('Flux (e$^{-}$/sec)')
if smoothtype == "b":
   exec("plt.title('%d: Raw with Boxcar Fit')" % kic)
elif smoothtype == "g":
   exec("plt.title('%d: Raw with Gaussian Fit')" % kic)
exec("plt.savefig('kic%d_raw.png')" % kic)

# concatenation for smoothed data

time = np.zeros(0)
sap_flux = np.zeros(0)

for i, q in enumerate(quarters):
   exec("time = np.append(time, time_%d)" % q)
   exec("sap_flux = np.append(sap_flux, sap_flux2_%d)" % q)

# scan for outliers
print ' '
while True:
   try:
      inp = int(raw_input('--> Clipping level: '))
      break
   except ValueError:
      print '--> Please enter an integer'
print ' '
clip = inp * np.std(sap_flux)
meanflux = np.mean(sap_flux)

upperbound = meanflux + clip
lowerbound = meanflux - clip

colours = np.zeros(sap_flux.size)

for i, flux in enumerate(sap_flux):
   if flux < upperbound and flux > lowerbound:
      colours[i] = 1

clipped_flux = []
clipped_time = []

# smoothed data
plt.figure(2)
plt.plot(time, sap_flux, 'bo', markersize=3)
plt.xlabel('Time (d)')
plt.ylabel('Fractional Intensity')
exec("plt.title('%d')" % kic)
for i, colour in enumerate(colours):
   if colour == 1:
      clipped_flux.append(sap_flux[i])
      clipped_time.append(time[i])
plt.plot(clipped_time, clipped_flux, 'ro', markersize=3)
exec("plt.savefig('kic%d_smooth.png')" % kic)

# export smoothed and clipped data as .dat file
exportblend = np.array([clipped_time, clipped_flux])
exportblend = np.transpose(exportblend)
exec("np.savetxt('kic%d_lc.dat', exportblend, delimiter=' ', header='Smoothed and clipped light curve for KIC%d')" % (kic, kic))

show = raw_input('--> Show plots now? (y/n) ')
while show != "y" and show != "n":
   show = raw_input('--> Please enter y or n: ')
print ' '
if show == "y":
   plt.show()
else:
   pass
