# pychimera
Data processing and analysis pipeline for CHIMERA instrument

Python package to reduce, process and analyze science images from the CHIMERA instrument on Hale 200-inch telescope at Palomar.

* Author: Navtej Singh
* Contact: navtej@astro.caltech.edu
* Website: http://tauceti.caltech.edu/chimera
* Organization: Astronomy@Caltech www.astro.caltech.edu
* Description: Python package to process CHIMERA instrument images

* Following requirements should be met to run the code.
  + Numpy >= 1.8
  + pyfits >= 0.14
  + pyraf >= 2.1
  + photutils >= 0.2

* To install pychimera -
  + Clone the pychimera repository to local machine
  + Run setup.py script to install the package 
    * python setup.py install
  
* Example scripts using chimera python modules is included in script directory. Three scripts are included -
  + animation.py : Generate movie from the frames of CHIMERA image cube
  + reduce.py : De-bias and flat field the CHIMERA science images
  + photometry.py : Perform aperture photometry of a star in CHIMERA science images and generate a light curve
