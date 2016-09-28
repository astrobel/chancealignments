import numpy as np
import scipy as sp
from astropy.io import fits as pyfits
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import lomb
import smoothing 
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib as mpl
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

quarterlist = []

for i in os.listdir(os.getcwd()):
   if i.endswith("lpd-targ.fits"):
      exec("hdulist_temp = pyfits.open('%s')" % i)
      parameters = hdulist_temp[0].header
      q = parameters['QUARTER']
      quarterlist.append(q)
      hdulist_temp.close()
   else:
      continue

print ' '
repquarter = 18 # placeholder, since there is 0 but there is no 18
while repquarter not in quarterlist:
   while True:
      try:
         repquarter = int(raw_input('--> Which quarter? (0-17) '))
         break
      except ValueError:
         print '--> Please enter an integer'
   if repquarter in quarterlist:
      break
   else:
      print '--> No data for this quarter'
      continue
print ' '
stitch = raw_input('--> Stitch with related quarters if possible? (y/n) ')
while stitch != "y" and stitch != "n":
   stitch = raw_input('--> Please type either y or n: ')
print ' '

stitchcounter = 1
grouping = repquarter % 4

for i in os.listdir(os.getcwd()):
   if i.endswith("lpd-targ.fits"):
      exec("hdulist_temp = pyfits.open('%s')" % i)
      parameters = hdulist_temp[0].header
      q = parameters['QUARTER']
      if q == repquarter:
         exec("hdulist1 = pyfits.open('%s')" % i)
         parameters = hdulist1[0].header
         kic = parameters['KEPLERID']
         channel = parameters['CHANNEL']
         table = hdulist1[1].data
         flux1 = table['FLUX']
         time1 = table['TIME']
         hd1 = hdulist1[1].header
         ysize1 = hd1['NAXIS2']
         table2 = hdulist1[2].data
         hd2 = hdulist1[2].header
         x = hd2['NAXIS1']
         y = hd2['NAXIS2']
         xsize = x * y
         temp2d = np.zeros((x, y))
         hdulist1.close()
      else:
         continue
   else:
      continue

if (channel%2) == 0:
   eo = 0
else:
   eo = 1

if stitch == "y":
   for i in os.listdir(os.getcwd()):
      if i.endswith("lpd-targ.fits"):
         exec("hdulist_temp = pyfits.open('%s')" % i)
         parameters = hdulist_temp[0].header
         q = parameters['QUARTER']
         parameters = hdulist_temp[2].header
         xtest = parameters['NAXIS1']
         ytest = parameters['NAXIS2']
         if (q % 4) == grouping and q != repquarter and q != 1 and xtest == x and ytest == y:
            stitchcounter += 1
            print 'yeehaw'
            exec("hdulist = pyfits.open('%s')" % i)
            parameters = hdulist[0].header
            hd = hdulist[1].header
            exec("ysize%d = hd['NAXIS2']" % stitchcounter)
            table = hdulist[1].data
            exec("flux%d = table['FLUX']" % stitchcounter)
            exec("time%d = table['TIME']" % stitchcounter)
            hdulist.close()
         else:
            continue
      else:
         continue
else:
   pass


### INPUT PARAMETERS ###

print ' '
print '---> Light curve smoothing <---'
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

print ' '
while True:
   try:
      inp = int(raw_input('--> Clipping level: '))
      break
   except ValueError:
      print '--> Please enter an integer'
print ' '
print ' '
print '--> Fourier transform <--'
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
print ' '


### ETC ###

# dynamic variable names
for (j, k), img in np.ndenumerate(temp2d):
   index = (k + 1) * j + (x - j) * k
   exec("pixel%d_flux = np.array(None)" % index)
   exec("pixel%d_time = np.array(None)" % index)

# filling the flux array
for l in range(stitchcounter):
   m = l + 1
   exec("second_flux%d = np.zeros([xsize, ysize%d])" % (m, m))
   exec("flux = flux%d" % m)
   for (i, j, k), val in np.ndenumerate(flux):
      index = (j + 1) * k + (x - k) * j
      exec("second_flux%d[index, i] = val" % m)

exp = raw_input('--> Export data for each pixel? (y/n) ')
while exp != "y" and exp != "n":
   exp = raw_input('--> Please enter y or n: ')
print ' '

for (j, k), img in np.ndenumerate(table2):
   index = (j + 1) * k + (x - k) * j
   if img == 0:
      pass
   else:

      flux_final = []
      time_final = []

      for l in range(stitchcounter):
         m = l + 1
         exec("flux_intermediate = second_flux%d[index,:]" % m)

         # creating a blend array to remove NaNs
         exec("blend = np.array([time%d, flux_intermediate])" % m)
         blend = np.transpose(blend)
         blend2 = np.ma.compress_rows(np.ma.fix_invalid(blend))

         exec("time_final%d = blend2[:,0]" % m)
         flux_intermediate = blend2[:,1]

         if smoothtype == "b": # boxcar smoothing
            exec("flux_final%d, smth_flux%d = smoothing.boxsmooth(time_final%d, flux_intermediate, kern)" % (m, m, m))
         elif smoothtype == "g": # gaussian smoothing
            exec("flux_final%d, smth_flux%d = smoothing.gausssmooth(time_final%d, flux_intermediate, kern)" % (m, m, m))

         exec("flux_final.extend(np.ndarray.tolist(flux_final%d))" % m)
         exec("time_final.extend(np.ndarray.tolist(time_final%d))" % m)

      flux_final = np.asarray(flux_final)
      time_final = np.asarray(time_final)

      exec("pixel%d_flux = flux_final" % index)
      exec("pixel%d_time = time_final" % index)

      exec("tempflux = pixel%d_flux" % index)
      exec("temptime = pixel%d_time" % index)

      clip = inp * np.std(tempflux)
      meanflux = np.mean(tempflux)
 
      upperbound = meanflux + clip
      lowerbound = meanflux - clip

      colours = np.zeros(tempflux.size)

      for i, flux in enumerate(tempflux):
         if flux < upperbound and flux > lowerbound:
            colours[i] = 1

      clipped_flux = []
      clipped_time = []
      for i, colour in enumerate(colours):
         if colour == 1:
            clipped_flux.append(tempflux[i])
            clipped_time.append(temptime[i])

      exec("pixel%d_flux = clipped_flux" % index)
      exec("pixel%d_time = clipped_time" % index)
 
      # export smoothed and clipped data as .dat file
      if exp == "y":
         exportblend = np.array([clipped_time, clipped_flux])
         exportblend = np.transpose(exportblend)
         exec("np.savetxt('kic%d_pixel%d_lc.dat', exportblend, delimiter=' ', header='Smoothed and clipped light curve for KIC%d TPF (please note that pixels are numbered on the rectanglar aperture from left to right, top to bottom, and inclusive of pixels with no data)')" % (kic, index+1, kic))
      else:
         pass

      # fourier transform
      frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
      if frequencies.shape != power_spectrum.shape:
         power_spectrum = np.concatenate((power_spectrum, [0]))
      power_spectrum = np.sqrt(power_spectrum)
      power_spectrum = 1e6 * ((power_spectrum / (np.median(power_spectrum))) - 1.) # to ppm

      exec("pixel%d_freq = frequencies" % index)
      exec("pixel%d_ps = power_spectrum" % index)


### PLOTTING ###

# light curves
fig = plt.figure(1)
gs = gridspec.GridSpec(y, x, wspace=0, hspace=0)
exec("plt.title('%d')" % kic)
plt.xlabel('Time (d)')
plt.ylabel('Fractional Intensity')

for (j, k), img in np.ndenumerate(table2):
   index = (j + 1) * k + (x - k) * j
   if img == 0:
      if eo == 0:
         ax = fig.add_subplot(gs[y - j - 1, x - k - 1])
      elif eo == 1:
         ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
   else:
      exec("flux = pixel%d_flux" % index)
      exec("time = pixel%d_time" % index)
      if eo == 0:
         ax = fig.add_subplot(gs[y - j - 1, x - k - 1])
      elif eo == 1:
         ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
      if img == np.amax(table2):
         plt.plot(time, flux, 'r-')
         plt.ylim(min(flux), max(flux)) #ymin=0)
         plt.xlim(min(time), max(time))
      else:
         plt.plot(time, flux, 'k-')
         plt.ylim(min(flux), max(flux)) #ymin=0)
         plt.xlim(min(time), max(time))

exec("plt.savefig('kic%d_pixelslc.png')" % kic)

# power spectra
fig = plt.figure(2)
gs = gridspec.GridSpec(y, x, wspace=0, hspace=0)
exec("plt.title('%d')" % kic)
plt.xlabel('Frequency ($\mu$Hz)')
plt.ylabel('Amplitude (ppm)')

for (j, k), img in np.ndenumerate(table2):
   index = (j + 1) * k + (x - k) * j
   if img == 0:
      if eo == 0:
         ax = fig.add_subplot(gs[y - j - 1, x - k - 1])
      elif eo == 1:
         ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
   else:
      exec("freq = pixel%d_freq" % index)
      exec("ps = pixel%d_ps" % index)
      if eo == 0:
         ax = fig.add_subplot(gs[y - j - 1, x - k - 1])
      elif eo == 1:
         ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
      if img == np.amax(table2):
         plt.plot(freq, ps, 'r-')
         plt.ylim(0, max(ps)) #ymin=0)
         plt.xlim(0, max(freq))
      else:
         plt.plot(freq, ps, 'k-')
         plt.ylim(0, max(ps)) #ymin=0)
         plt.xlim(0, max(freq))

exec("plt.savefig('kic%d_pixels.png')" % kic)

show = raw_input('--> Show plots now? (y/n) ')
while show != "y" and show != "n":
   show = raw_input('--> Please enter y or n: ')
print ' '
if show == "y":
   plt.show()
else:
   pass
