import numpy as np
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel


def boxsmooth(time, flux, kern):
   """Takes a time series and boxcar smoothes it,
         removing edge effects.
      time: ndarray, 1D
      flux: ndarray, 1D
      kern: boxcar smoothing
   """

   box_flux = convolve(flux, Box1DKernel(kern), boundary='extend')

   flux2 = flux / box_flux

   return flux2, box_flux


def gausssmooth(time, flux, kern):
   """Takes a time series and applies a Gaussian filter,
         removing edge effects.
      time: ndarray, 1D
      flux: ndarray, 1D
      kern: gaussian smoothing
   """

   gauss_flux = convolve(flux, Gaussian1DKernel(kern), boundary='extend')

   flux2 = flux / gauss_flux

   return flux2, gauss_flux # correct, fit


def polysmooth(time, flux, degree):
   """Takes a time series and fits a polynomial
      of input degree, then divides by polynomial.
   """

   pfit = np.polyfit(time, flux, degree)
   fit = np.poly1d(pfit)
   xfit = np.linspace(min(time), max(time), len(time))

   flux2 = flux / fit(xfit)

   return flux2, fit, xfit
