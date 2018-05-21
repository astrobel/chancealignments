import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits as pyfits
from astropy import wcs
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
from astropy.utils.exceptions import AstropyWarning
import lomb
import smoothing
import translate as tr
import quarters as qs
import nancleaner as nc
import os, sys, warnings

warnings.simplefilter('ignore', category=AstropyWarning)

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
      r'\usepackage{helvet}',
      r'\usepackage[EULERGREEK]{sansmath}',
      r'\sansmath'
]
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True

def fft(time, flux, ofac, hifac):
   freq,power, nout, jmax, prob = lomb.fasper(time, flux, ofac, hifac)
   convfactor = (1. / (60 * 60 * 24)) * (10 ** 6)
   uHzfreq = freq * convfactor #11.57, conversion c/d to mHz
   return uHzfreq, power

# pick a quarter
repquarter = qs.getquarter()

# grab the quarter
for i in os.listdir(os.getcwd()):
   if i.endswith("lpd-targ.fits"):
      exec("hdulist_temp = pyfits.open('%s')" % i)
      things = hdulist_temp[0].header
      q = things['QUARTER']
      if q == repquarter:
         exec("hdulist1 = pyfits.open('%s')" % i)
         parameters = hdulist1[0].header
         kic = parameters['KEPLERID']
         channel = parameters['CHANNEL']
         table = hdulist1[1].data
         flux1 = table['FLUX']
         time1 = table['TIME']
         hd1 = hdulist1[1].header
         w = wcs.WCS(hd1, keysel=['binary'])
         hd2 = hdulist1[2].header
         x = hd2['NAXIS1']
         y = hd2['NAXIS2']
         refx = hd2['CRPIX1']
         refy = hd2['CRPIX2']
         refra = hd2['CRVAL1']
         refdec = hd2['CRVAL2']
         hdulist1.close()
      else:
         continue
   else:
      continue

if (channel%2) == 0:
   eo = 0
else:
   eo = 1

dim = len(flux1)
nans = []

# fluxnew = np.zeros((dim, y, x))

for i in range(dim):
   # fluxnew[i,:] = flux1[i,:]
   if np.isnan(flux1[i,:]).all() == True:
      # print 'ye'
      nans.append(i)

fluxnew = np.delete(flux1, nans, axis=0)

avgflux = np.mean(fluxnew, axis=0)

# read in light curve
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

binnum = 1000
binsize = foldper / binnum

phasedtimearray = np.zeros(len(times))
finalampls = np.zeros(binnum)
amplcounts = np.zeros(binnum)

for i, val in enumerate(times):
   phasedtime = val % foldper
   newphasedtime = tr.translate(phasedtime, 0, foldper, 0, 1)
   phasedtimearray[i] = newphasedtime
   bindex = int((phasedtime - (phasedtime % binsize)) / binsize - 1)
   finalampls[bindex] += ampls[i]
   amplcounts[bindex] += 1
   
finalampls = np.divide(finalampls, amplcounts)

finaltimes = np.histogram(phasedtimearray, bins=binnum-1, range=(0,1))

nanlist = []
for i, val in enumerate(finalampls):
   if np.isnan(val) == True:
      nanlist.append(i)

phasetime = np.delete(finaltimes[1], nanlist)
phaseflux = np.delete(finalampls, nanlist)

maxpos = np.argmax(phaseflux) / np.float64(len(phasetime))
minpos = np.argmin(phaseflux) / np.float64(len(phasetime))

# folding the pixel fluxes
time1 = time1 % foldper
for i, val in enumerate(time1):
   time1[i] = tr.translate(val, 0, foldper, 0, 1)

# find the fluxes for differencing
timeflags = np.zeros(len(time1))
tolerance = 0.05

for i, time in enumerate(time1):
   if minpos - tolerance <= 0:
      if time < minpos + tolerance or time > (minpos - tolerance)%1:
         timeflags[i] = -1
   elif maxpos + tolerance >= 1:
      if time < (minpos + tolerance)%1 or time > minpos - tolerance:
         timeflags[i] = -1
   else:
      if time < minpos + tolerance and time > minpos - tolerance:
         timeflags[i] = -1

   if maxpos - tolerance <= 0:
      if time < maxpos + tolerance or time > (maxpos - tolerance)%1:
         timeflags[i] = 1
   elif maxpos + tolerance >= 1:
      if time < (maxpos + tolerance)%1 or time > maxpos - tolerance:
         timeflags[i] = 1
   else:
      if time < maxpos + tolerance and time > maxpos - tolerance:
         timeflags[i] = 1

# establish number of frames to be used around the maxima and minima
highflags = np.where(timeflags>0)[0]
lowflags = np.where(timeflags<0)[0]
highdim = len(highflags)
lowdim = len(lowflags)

# 3d flux arrays to fill with frames
fluxhigh = np.zeros((highdim, y, x))
fluxlow = np.zeros((lowdim, y, x))

# get rid of nans
highnans = []
lownans = []
for i in range(highdim):
   fluxhigh[i,:] = flux1[highflags[i]]
   if np.isnan(fluxhigh[i,:]).all() == True:
      highnans.append(i)
for i in range(lowdim):
   fluxlow[i,:] = flux1[lowflags[i]]
   if np.isnan(fluxlow[i,:]).all() == True:
      lownans.append(i)

fluxhigh1 = np.delete(fluxhigh, highnans, axis=0)
fluxlow1 = np.delete(fluxlow, lownans, axis=0)

fluxdiff = np.abs(np.average(fluxhigh1, axis=0) - np.average(fluxlow1, axis=0))

# print fluxdiff

imgflux = np.flipud(fluxdiff)
if eo == 0:
   imgflux = np.fliplr(imgflux)
avgflux = np.flipud(avgflux)
if eo == 0:
   avgflux = np.fliplr(avgflux)


### PLOTTING ###

plt.figure(1)

fig, (avgimg, diffimg) = plt.subplots(1, 2) 

# exec("fig.suptitle('%d q%d', fontsize=20)" % (kic, repquarter))

crval = w.wcs.crval
north = crval + np.array([0, 6/3600.])
east = crval + np.array([ 6/3600., 0])

ncoords = np.vstack([crval, north])
ecoords = np.vstack([crval, east])
npixels = w.wcs_world2pix(ncoords , 0)
epixels = w.wcs_world2pix(ecoords , 0)
npixels[1, 1] = npixels[0, 1] - (npixels[1, 1] - npixels[0, 1]) # flip ud
epixels[1, 1] = epixels[0, 1] - (epixels[1, 1] - epixels[0, 1])
if eo == 0:
   npixels[1, 0] = npixels[0, 0] - (npixels[1, 0] - npixels[0, 0]) # flip lr
   epixels[1, 0] = epixels[0, 0] - (epixels[1, 0] - epixels[0, 0])
diffimg.plot(npixels[:,0], npixels[:,1], color='#0cb5ed')
diffimg.plot(epixels[:,0], epixels[:,1], '--', color='#0cb5ed')

left = avgimg.imshow(avgflux, cmap='YlOrRd')
avgimg.set_title('Average')
left.set_interpolation('nearest')
avgimg.set_xlim(-0.5, x-0.5)
avgimg.set_ylim(y-0.5, -0.5)
avgimg.plot(npixels[:,0], npixels[:,1], color='#0cb5ed')
avgimg.plot(epixels[:,0], epixels[:,1], '--', color='#0cb5ed')
cbar = fig.colorbar(left, ax=avgimg)

left.axes.get_xaxis().set_ticklabels([])
left.axes.get_yaxis().set_ticklabels([])
left.axes.get_xaxis().set_ticks([])
left.axes.get_yaxis().set_ticks([])

avgimg.plot([25, 25], [25, 55], '-', color='#0cb5ed')
avgimg.plot([25, 55], [25, 25], '--', color='#0cb5ed')

right = diffimg.imshow(imgflux, cmap='YlOrRd')
diffimg.set_title('Difference')
right.set_interpolation('nearest')
diffimg.set_xlim(-0.5, x-0.5)
diffimg.set_ylim(y-0.5, -0.5)
cbar2 = fig.colorbar(right, ax=diffimg)

right.axes.get_xaxis().set_ticklabels([])
right.axes.get_yaxis().set_ticklabels([])
right.axes.get_xaxis().set_ticks([])
right.axes.get_yaxis().set_ticks([])

if eo == 1:
   diffimg.plot((x - refx) - 0.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)
   avgimg.plot((x - refx) - 0.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)
elif eo == 0:
   diffimg.plot(refx - 1.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)
   avgimg.plot(refx - 1.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)

plt.tight_layout()
fig.set_size_inches(7.5, 4.5)
exec("plt.savefig('kic%dq%dimgs.png')" % (kic, repquarter))

print ' '
show = raw_input('--> Show plot now? (y/n) ')
while show != "y" and show != "n":
   show = raw_input('--> Please enter y or n: ')
print ' '
if show == "y":
   plt.show()
else:
   pass
