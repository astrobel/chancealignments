import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits as pyfits
from astropy import wcs
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
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True

def fft(time, flux, ofac, hifac):
	freq,power, nout, jmax, prob = lomb.fasper(time, flux, ofac, hifac)
	convfactor = (1. / (60 * 60 * 24)) * (10 ** 6)
	uHzfreq = freq * convfactor #11.57, conversion c/d to mHz
	return uHzfreq, power

# pick a quarter
quarterlist = []

for i in os.listdir(os.getcwd()):
	if i.endswith("lpd-targ.fits"):
		exec("hdulist_temp = pyfits.open('%s')" % i)
		things1 = hdulist_temp[0].header
		q1 = things1['QUARTER']
		quarterlist.append(q1)
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

for i in range(dim):
	if np.isnan(flux1[i,:]).all() == True:
		nans.append(i)

fluxnew = np.delete(flux1, nans, axis=0)

avgflux = np.mean(fluxnew, axis=0)

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
   print ' '
   frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
   power_spectrum = np.concatenate((power_spectrum, [0]))
   power_spectrum = power_spectrum * 4 * np.var(clipped_flux) / clipped_flux.size
   power_spectrum = np.sqrt(power_spectrum)
   power_spectrum *= 1e6 #* ((power_spectrum / (np.median(power_spectrum))) - 1.) # to ppm

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

# binning
while True:
	try:
		binnum = int(raw_input('--> Number of bins: (integer) '))
		break
	except ValueError:
		print '--> Please enter an integer'
binnum = np.float64(binnum)
print ' '

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

nanlist = []
for i, val in enumerate(flux_sum):
	if np.isnan(val) == True:
		nanlist.append(i)

phasetime = np.delete(time_binned, nanlist)
phaseflux = np.delete(flux_sum, nanlist)

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
	if time < maxpos + tolerance and time > maxpos - tolerance:
		timeflags[i] = 1
	elif time < minpos + tolerance and time > minpos - tolerance:
		timeflags[i] = -1

highflags = np.where(timeflags>0)[0]
lowflags = np.where(timeflags<0)[0]
highdim = len(highflags)
lowdim = len(lowflags)
highnans = []
lownans = []

fluxhigh = np.zeros((highdim, y, x))
fluxlow = np.zeros((lowdim, y, x))

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

imgflux = np.flipud(fluxdiff)
if eo == 0:
	imgflux = np.fliplr(imgflux)
avgflux = np.flipud(avgflux)
if eo == 0:
	avgflux = np.fliplr(avgflux)


### PLOTTING ###

plt.figure(1)

fig, (ukirt, kepler) = plt.subplots(1, 2) 

left = kepler.imshow(imgflux, cmap='YlOrRd')
kepler.set_title('Difference')
left.set_interpolation('nearest')
kepler.set_xlim(-0.5, x-0.5)
kepler.set_ylim(y-0.5, -0.5)

left.axes.get_xaxis().set_ticklabels([])
left.axes.get_yaxis().set_ticklabels([])
left.axes.get_xaxis().set_ticks([])
left.axes.get_yaxis().set_ticks([])

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
kepler.plot(npixels[:,0], npixels[:,1], color='#0cb5ed')
kepler.plot(epixels[:,0], epixels[:,1], '--', color='#0cb5ed')

if eo == 1:
   kepler.plot((x - refx) - 0.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)
   ukirt.plot((x - refx) - 0.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)
elif eo == 0:
   kepler.plot(refx - 1.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)
   ukirt.plot(refx - 1.5, (y - refy) - 0.5, '*', color='#0cb5ed', ms=10)

right = ukirt.imshow(avgflux, cmap='YlOrRd')
ukirt.set_title('Average')
right.set_interpolation('nearest')
ukirt.set_xlim(-0.5, x-0.5)
ukirt.set_ylim(y-0.5, -0.5)
ukirt.plot(npixels[:,0], npixels[:,1], color='#0cb5ed')
ukirt.plot(epixels[:,0], epixels[:,1], '--', color='#0cb5ed')

right.axes.get_xaxis().set_ticklabels([])
right.axes.get_yaxis().set_ticklabels([])
right.axes.get_xaxis().set_ticks([])
right.axes.get_yaxis().set_ticks([])

ukirt.plot([25, 25], [25, 55], '-', color='#0cb5ed')
ukirt.plot([25, 55], [25, 25], '--', color='#0cb5ed')

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
