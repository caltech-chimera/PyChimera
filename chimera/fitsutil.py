from __future__ import division


import pyfits

import chimera


def fitsread(imgname, header = False):
    """
    Read CHIMERA FITS image cube.
    
    Parameters
    ----------
    image : string
        FITS image name
        
    header : boolean
        Return FITS image header?
        
    Returns
    -------
    img_data : numpy array
        2D or 3D numpy array
    """
    try:
        if header:
            img_data, header = pyfits.getdata(imgname, ignore_missing_end = True, header = True)
            return img_data, header
        else:
            img_data = pyfits.getdata(imgname, ignore_missing_end = True)
            return img_data
    except IOError:
        print "FITSREAD: Unable to open FITS image %s" %imgname
    
    return
    
    
def fitswrite(img, imgname, **kwargs):
    """
    Write FITS image to disk.
    
    Parameters
    ----------
    img : numpy array
        2D or 3D numpy array
        
    imgname : string
        Name of the output FITS image
        
    Optional Keywords
    -----------------
    header : pyFITS header
        FITS header object
        
    Return
    ------
    None
    """
    try:
        if kwargs.has_key('header'):
            hdu = pyfits.PrimaryHDU(img, header = kwargs['header'])
        else:
            hdu = pyfits.PrimaryHDU(img)
        hdu.writeto(imgname)
    except IOError:
        print "FITSWRITE: Unable to write FITS image %s. Stopping." %imgname
    
    return
    
    
def fitshead(imgname):
    """
    Read CHIMERA FITS image header.
    
    Parameters
    ----------
    image : string
        FITS image name
        
    Returns
    -------
    img_header : python dictionary
        Dictionary of image header keywords
    """
    try:
        img_header = pyfits.getheader(imgname, ignore_missing_end = True)
        return img_header
    except IOError:
        print "FITSHEAD: Unable to open FITS image %s. Stopping." %imgname
    
    return
    
    
def fitscombine(imgname, combine = "average", nframes = 100, outfile = ""):
    """
    Combine FITS image frames of CHIMERA instrument 3D image cubes.
    
    Parameters
    ----------
    imgname : string
        FITS image cube name
        
    combine : string
        Combine type - average, sum or median
    
    nframes : int
        Number of frames to combine
        
    outfile : string
        Name of the output FITS image
        
    Returns
    -------
    None
    """
    img_data, img_header = fitsread(imgname, header = True)
    
    combine_image = chimera.imcombine(img_data, combine, nframes)
    
    if outfile == "":
        fitswrite(combine_image, imgname.replace(".fits", "_comb.fits"), header = img_header)
    else:
        fitswrite(combine_image, outfile, header = img_header)
        
    return
