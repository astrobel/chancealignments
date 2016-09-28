--> An automated pipeline to detect chance alignments in Kepler light curves and target pixel files <--
--> Isabel C., 2016 <--
--> github.com/astrobel <--

Features:
--> Fully automated; no need to modify files
--> Lots of user input parameters to customise analysis
--> Error handling for all user input
--> Inbuilt methods for dealing with the odd and even Kepler CCD readout problem
--> Option to show all plots at the end of running each script
--> Code to visually compare Kepler pixel apertures to UKIRT (UK Infrared Telescope) images
--> !! NEW !! Difference imaging code to locate contaminant source and false positives
--> !! NEW !! Options for phasing light curves on peaks that do not have the highest amplitude in Fourier space

You will need:
--> All your Kepler llc.fits and lpd-targ.fits files in one folder
--> 1_smoothing.py
--> 2_spectrum.py
--> 3_phase.py
--> 4_pixels.py
--> 5_comparison.py
-----> A UKIRT FITS image of the target: http://wsa.roe.ac.uk:8080/wsa/getImage_form.jsp (requires login as outlined on http://keplergo.arc.nasa.gov/ToolsUKIRT.shtml)
--> 6_locations.py
-----> A UKIRT FITS image of the target (as above)
-----> kcoords.dat: a file of all Kepler objects within 1' of the target: http://archive.stsci.edu/kepler/kepler_fov/search.php
-----> ucoords.dat: a file of all UKIRT objects within 1' of the target: http://wsa.roe.ac.uk:8080/wsa/region_form.jsp (requires login)
--> 7_difference.py
--> lomb.py
--> smoothing.py
--> translate.py
--> Absolutely no other FITS files!

How to set up the files for part 6:
--> kcoords.dat
-----> Column 1: RA
-----> Column 2: DEC
-----> Column 3: Kepler magnitude
--> ucoords.dat
-----> Column 1: RA
-----> Column 2: DEC
-----> Column 3: J magnitude

PART 1: Preparing the light curves
--> Run 1_smoothing.py
--> Prompt: Boxcar or Gaussian smoothing? (b/g)
-----> (Explain different types of smoothing)
--> Prompt: Smoothing kernel: (integer)
-----> This is the same for boxcar or Gaussian smoothing. I use 100 width with Gaussian smoothing for best results.
--> Prompt: Clipping level: (integer)
-----> This is to remove outliers from the final light curve. Outliers will be clipped above an inputted number of standard deviations from the mean. I find 3 sigma is sufficient.
--> Prompt: Show plots now? (y/n)
-----> Matplotlib has a built in option plot.show() which displays interactive plots. This way, it's easy to zoom in and hover over a data point to measure it. Either way, plots will be saved as static PNG files.

PART 2: Obtaining an amplitude spectrum
--> Run 2_spectrum.py
--> Prompt: Oversampling factor: (integer)
-----> For this, a lower oversampling factor is sufficient. I use 5.
--> Prompt: Nyquist range factor: (float)
-----> In general, this will be 1, unless you want to sample beyond the Nyquist frequency, which this code is not particularly designed to do. This can also be set to a float less than 1 to enable closer inspection of lower frequencies.
--> Prompt: Show plot now? (y/n)

PART 3: Phasing the light curve
--> Run 3_phase.py
--> Prompt: Use highest peak for folding? (y/n)
-----> If y, this starts by taking another Fourier transform so it can find the location of the highest peak for folding.
--------> Prompt: Oversampling factor: (integer)
-----------> For this, higher detail is optimal. I use an oversampling factor of 20.
--------> Prompt: Nyquist range factor: (float)
--------> Prompt: Highest frequency to plot in microHertz: (float)
--------> This code produces two plots, the first of which is an optionally zoomed-in amplitude spectrum to show more detail. It will plot up to whichever frequency you enter here, which can be chosen based on the spectrum obtained in part 2.
-----> If n:
--------> Prompt: Enter folding frequency: (microHertz)
-----------> Can be obtained manually from 2_spectrum.py
--> Prompt: Smooth or bin the phase curve? (s/b)
-----> Smoothing: As in part 1, choose between boxcar and Gaussian smoothing
-----> Binning: Prompt: Number of bins (integer)
--> Prompt: Show plots now? (y/n)
-----> The second plot shows the phase curve, with red points indicating the inputted number of bins, and blue points indicating one tenth that number of bins.

PART 4: Analysing individual pixels
--> Run 4_pixels.py
--> Prompt: Which quarter? (0-17)
-----> In many cases, Kepler targets won't have all quarters available. If you enter a number of a quarter that doesn't exist, the code will inform you of this and prompt you for another.
--> Prompt: Stitch with related quarters if possible? (y/n)
-----> To get greater detail from pixel data, it is useful to combine quarters that fall on the same module and have apertures of the same size. The code searches for these instances and uses all available quarters. In many cases, there won't be any that match.
--> Prompt: Boxcar or Gaussian smoothing (b/g)
--> Prompt: Smoothing kernel: (integer)
--> Prompt: Clipping level: (integer)
--> Prompt: Oversampling factor: (integer)
-----> This doesn't need as much detail as in part 3, so I return to oversampling by a factor of 5.
--> Prompt: Nyquist range factor: (float)
--> Prompt: Export data for each pixel? (y/n)
-----> This exports the light curve for each individual pixel, which can be used for more detailed Fourier analysis. Information about how to identify the pixels is contained in the header of each exported light curve.
--> Prompt: Show plots now? (y/n)
-----> Plots produced show the pixel aperture with light curves and amplitude spectra plotted in each pixel.

PART 5: Comparing the pixel aperture to a UKIRT image
--> Run 5_comparison.py
--> Prompt: Which quarter? (0-17)
-----> For the best comparison, use the same as in part 4.
--> Prompt: Plot the reference pixel on the Kepler image? (y/n)
-----> The reference pixel refers to the location on the Kepler aperture that corresponds to the provided RA and DEC, which will be printed in the the terminal if the user enters y.
--> Prompt: Show plot now? (y/n)
-----> The Kepler aperture is plotted in the left panel, and the UKIRT image in the right. In both cases a compass rose is shown with North as a solid line and East as a dashed like. The North arm is 6" to scale in each image.

PART 6: Looking for possible contaminating targets within 1' of the target
--> Run 6_locations.py
--> Prompt: Show plot now? (y/n)
-----> Two plots are produced, one which shows the UKIRT image with the target star circled, and one which has the same base but includes an overlay of all nearby targets in the UKIRT (circles) and Kepler (stars) catalogues. J magnitudes are indicated by a colour scale.

PART 7: Difference imaging
--> Run 7_difference.py
--> Prompt: Which quarter? (1-17)
--> For now, the code is the same as 3_phase.py, as it requires you to choose (manually or otherwise) a peak for phase folding. The code takes light curve data up to 10% either side of each maximum and minimum in the phase curve and subtracts the average of the flux around the minima from the average flux around the maxima to create the difference image.
--> Prompt: Number of bins: (integer)
-----> It is best to choose this by inspection of the phase curve produced in 3_phase.py, in order to get the most accurate highest and lowest phased flux to choose maxima and minima.
--> Prompt: Show plot now? (y/n)
