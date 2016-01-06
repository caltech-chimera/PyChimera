#!/usr/bin/env python

"""
    --------------------------------------------------------------------------
    Routine to generate animation of CHIMERA image cubes.
    
    Usage: python animation.py [options] image_cube
        
                                
    Author:
        Navtej Saini

    Organization:
        Caltech, Pasadena, CA, USA

    Version:
        23 December 2015     0.1dev     Initial implementation 
    --------------------------------------------------------------------------        
"""


import chimera
from optparse import OptionParser


# Matplotlib colormaps (Takne from http://matplotlib.org/examples/color/colormaps_reference.html)
cmaps = ['Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
         'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 
         'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'bone', 'cool', 'copper',
         'gist_heat', 'gray', 'hot', 'pink', 'spring', 'summer', 'winter',
         'BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 
         'RdYlBu', 'RdYlGn', 'Spectral', 'seismic', 'Accent', 'Dark2', 'Paired', 
         'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3', 'gist_earth', 'terrain', 
         'ocean', 'gist_stern', 'brg', 'CMRmap', 'cubehelix', 'gnuplot', 
         'gnuplot2', 'gist_ncar', 'nipy_spectral', 'jet', 'rainbow',
         'gist_rainbow', 'hsv', 'flag', 'prism']
         

if __name__ == "__main__":
    usage = "Usage: python %prog [options] image"
    description = "Description. Utility to generate animation of multi-extension CHIMERA FITS images."
    parser = OptionParser(usage = usage, version = "%prog 0.1dev", description = description)
    parser.add_option("-n", "--nframes", dest = "nframes", metavar="NFRAMES",
                    action="store", help = "Number of frames (default is 200)",
                    default = 200
    )
    parser.add_option("-y", "--vmin", dest = "vmin",
                    action='store', metavar="VMIN", help = "Minimum pixel value (default is 400)",
                    default = 400
    )
    parser.add_option("-x", "--vmax", dest = "vmax",
                    action='store', metavar="VMAX", help = "Maximum pixel value (default is 1000)",
                    default = 1500
    )
    parser.add_option("-s", "--scale", dest = "scale",
                    action='store', metavar="SCALE", help = "Image scaling (default is linear)",
                    choices = ['linear', 'log'], default = "linear"
    )     
    parser.add_option("-c", "--cmap", dest = "cmap",
                    action='store', metavar="CMAP", help = "Animation colormap (default is jet)",
                    choices = cmaps, default = "jet"
    ) 
    parser.add_option("-f", "--fps", dest = "fps",
                    action='store', metavar="FPS", help = "Animation frames per second (default is 30)",
                    default = 30
    )   
    parser.add_option("-o", "--outtype", dest = "outtype",
                    action='store', metavar="OUTTYPE", help = "Format of output animation (default is mp4)",
                    choices = ["mp4", "gif"], default = "mp4"
    )
    
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("Incorrect number of arguments")

    fname = chimera.animate(args[0], nframes = int(options.nframes), vmin = float(options.vmin), vmax = float(options.vmax), scale = options.scale, cmap = options.cmap, fps = options.fps, outtype = options.outtype)
    
    print "  Animation %s generated." %(fname)