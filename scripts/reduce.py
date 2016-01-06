#!/usr/bin/env python

"""
    --------------------------------------------------------------------------
    Routine to reduce raw science images from CHIMERA instrument.
    
    Usage: python reduce.py [options] image
        
                                
    Author:
        Navtej Saini

    Organization:
        Caltech, Pasadena, CA, USA

    Version:
        20 December 2015     0.1dev     Initial implementation 
    --------------------------------------------------------------------------        
"""

import os, sys
from StringIO import StringIO
from optparse import OptionParser

import chimera


def process(sci_files, bias_file, flat_file, nskip, threshold):
    """
    Entry point function to process science images.
    
    Parameters
    ----------
    sci_files : string
        Science image file names
        
    bias_file : string
        Bias image file name
        
    flat_file : string
        Flat field file name
        
    skip : int
        Number of bias frames to skip before averaging the frames. Default is 0.
        
    threshold : float
        Threshold for normalized fat field (value between 0 and 1.0). 
        Default is 0.8.
        
    Returns
    -------
    None 
    """
    print "REDUCE: CHIMERA Image Reduction Rotuine"
    
    nskip = int(nskip)
    threshold = float(threshold)
         
    # Check if input is a string of FITS images or a text file with file names
    if sci_files[0] == "@":
        infile = sci_files[1:]
        
        if not os.path.exists(infile):
            print "  Not able to locate file %s" %infile
        
        image_cubes = []
        with open(infile, "r") as fd:
            for line in fd.readlines():
                if len(line) > 1:
                    image_cubes.append(line.replace("\n", ""))
    else:
        image_cubes = sci_files.split(",")

    # Read bias and flat field images
    bias_image = chimera.fitsread(bias_file)
    flat_image = chimera.fitsread(flat_file)

    # Create master bias image
    print "  Generating master bias image"
    master_bias_image = chimera.masterbias(bias_image)

    # Create normalized flat field
    print "  Generating normalized flat field image"
    master_flat_image = chimera.masterflat(flat_image, master_bias_image)
        
    ncubes = len(image_cubes)
    for i in range(ncubes):
        sci_file = image_cubes[i]
        
        print "  Reducing science image : ", sci_file

        sci_image, header = chimera.fitsread(sci_file, header = True)

        # Reduced the science frames
        sci_red_image, sci_avg_image = chimera.imreduce(sci_image, master_bias_image, master_flat_image)

        # Write the reduced and average FITS image
        red_file = sci_file.replace('.fits', '_final.fits')
        avg_file = sci_file.replace('.fits', '_avg.fits')
        
        if os.path.exists(red_file):
            os.remove(red_file)
        if os.path.exists(avg_file):
            os.remove(avg_file)
        
        chimera.fitswrite(sci_red_image, red_file, header = header)
        chimera.fitswrite(sci_avg_image, avg_file, header = header)

        print "  Reduced science image : ", red_file

    return


if __name__ == "__main__":
    usage = "Usage: python %prog [options] sci_image bias_image flat_image"
    description = "Description. Utility to reduce raw science CHIMERA instrument images."
    parser = OptionParser(usage = usage, version = "%prog 0.1dev", description = description)
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default = False,
                      help = "print result messages to stdout"
                      )
    parser.add_option("-q", "--quiet",
                    action="store_false", dest="verbose", default = True,
                    help = "don't print result messages to stdout"
                    )
    parser.add_option("-s", "--skip", dest = "skip",
                    action='store', metavar="SKIP", help = "Number of frames to skip for average master bias (default is 0)",
                    default = 0
                    )
    parser.add_option("-t", "--threshold", dest = "threshold",
                    action='store', metavar="THRESHOLD", help = "Threshold for normalized flatfields (default is 0.8)",
                    default = 0.8
                    )
                                        
    (options, args) = parser.parse_args()  
    if len(args) != 3:
        parser.error("REDUCE: Incorrect number of arguments")
        
    # Check verbosity
    if not options.verbose:
        output = StringIO()
        old_stdout = sys.stdout
        sys.stdout = output
 
    process(args[0], args[1], args[2], options.skip, options.threshold)    

    # Reset verbosity
    if not options.verbose:
        sys.stdout = old_stdout