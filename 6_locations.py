import numpy as np
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import translate as tr
import os

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
      r'\usepackage{helvet}',
      r'\usepackage[EULERGREEK]{sansmath}',
      r'\sansmath'
]

for i in os.listdir(os.getcwd()):
   if i.endswith("llc.fits"):
      exec("hdulist = pyfits.open('%s')" % i)
      parameters = hdulist[0].header
      kic = parameters['KEPLERID']
      centra = parameters['RA_OBJ']
      centdec = parameters['DEC_OBJ']
      pass
   else:
      continue


### UKIRT IMAGE ###

for i in os.listdir(os.getcwd()):
   if i.endswith(".fits"):
      if not i.endswith("llc.fits") and not i.endswith("lpd-targ.fits"):
         exec("hdulist = pyfits.open('%s')" % i)
         break
   else:
      continue

flux2 = hdulist[1].data

print ' '
adj = raw_input('--> Adjust the UKIRT image orientation? (y/n) ')
while adj != "y" and adj != "n":
   adj = raw_input('--> Please enter y or n: ')
if adj == "y":
   rotby = 1
   while (rotby % 90) != 0:
      while True:
         try:
            rotby = int(raw_input('--> Enter a multiple of 90 degrees for rotation: '))
            break
         except ValueError:
            print '--> Please enter an integer'
      if (rotby % 90) == 0:
         break
      else:
         continue
   rotby = np.float64(rotby)
   rotby /= 90
   flux2 = np.rot90(flux2, rotby)
else:
   pass
print ' '


### READING IN SOME OTHER FILES ###
uk = np.loadtxt('ucoords.dat')
kep = np.loadtxt('kcoords.dat')

kepra = kep[:,0]
kepdec = kep[:,1]
kepmag = kep[:,2]
ukra = uk[:,0]
ukdec = uk[:,1]
ukmag =  uk[:,2]

for i, val in enumerate(kepmag):
   j = ukmag[i]
   if val <= 16.7:
      kepmag[i] = j - 398.04666 + 149.08127*j - 21.952130*(j**2) + 1.5968619*(j**3) - 0.057478947*(j**4) + 0.00082033223*(j**5)
   elif val > 16.7:
      kepmag[i] = j + 0.1918 + 0.08156*j


### PLOTTING ###

plt.figure(1) 

ukirt = plt.imshow(flux2, cmap='pink', norm=LogNorm())
exec("plt.title('%d')" % kic)
ukirt.set_interpolation('bilinear')
plt.xlim(300, 0)
plt.ylim(0, 300)

ukirt.axes.get_xaxis().set_ticklabels([])
ukirt.axes.get_yaxis().set_ticklabels([])
ukirt.axes.get_xaxis().set_ticks([])
ukirt.axes.get_yaxis().set_ticks([])

plt.plot([25, 25], [25, 55], '-', color='#00ff8c')
plt.plot([25, 55], [25, 25], '-', color='#00ff8c')

plt.text(30, 60, 'N', color='#00ff8c')
plt.text(70, 20, 'E', color='#00ff8c')

main = plt.Circle((149, 149), 15, color='#00ff8c', fill=False)
plt.gca().add_artist(main)
# opt_ = plt.Circle((x, y), 8, color='#00ff8c', fill=False)
# plt.gca().add_artist(opt_)
# plt.text(x-12, y, 'KIC', color='#00ff8c')

exec("plt.savefig('kic%dloc.png')" % kic)

cutoffra = 1/100. + centdec/100000 # i suspect the image is not precisely 1 arcmin in RA???
cutoffdec = 1/120. # 1/2 arcmin in degrees

for i, val in enumerate(ukra):
    ukra[i] = tr.translate(val, centra-cutoffra, centra+cutoffra, 0, 300)
for i, val in enumerate(kepra):
    kepra[i] = tr.translate(val, centra-cutoffra, centra+cutoffra, 0, 300)
for i, val in enumerate(ukdec):
    ukdec[i] = tr.translate(val, centdec-cutoffdec, centdec+cutoffdec, 0, 300)
for i, val in enumerate(kepdec):
    kepdec[i] = tr.translate(val, centdec-cutoffdec, centdec+cutoffdec, 0, 300)

plt.scatter(ukra, ukdec, s=50, c=ukmag, cmap='gist_rainbow', linewidths=0, alpha = 0.5)
plt.scatter(kepra, kepdec, marker=u'*', s=60, c=kepmag, cmap='gist_rainbow', linewidths=0.5, edgecolor='white')#, alpha = 0.7)

cbar = plt.colorbar()
cbar.ax.invert_yaxis()

exec("plt.savefig('kic%dloc1.png')" % kic)

show = raw_input('--> Show plots now? (y/n) ')
while show != "y" and show != "n":
   show = raw_input('--> Please enter y or n: ')
print ' '
if show == "y":
   plt.show()
else:
   pass
