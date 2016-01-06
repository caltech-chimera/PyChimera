from __future__ import division


import numpy as np


def imfill(image, xmin, xmax, ymin, ymax, value = -100):
    """
    Fill image section with some constant value.
    
    Parameters
    ----------
    fname : numpy array
        Image could be a 2D image or a 3D image cube.
        
    xmin, xmax : int
        Pixel coordinate along x-axis
        
    ymin, ymax : int
        Pixel coordinates along y-axis
        
    value : int
        Value to fill the image section with.
        
    Returns
    -------
    outfile : string
        FITS image with image section filled.
    """
    if len(image.shape) == 2:
        image[ymin:ymax, xmin:xmax] = value
    elif len(image.shape) == 3:
        image[:,ymin:ymax,xmin:xmax] = value
    
    return
    
    
def imshow(image, ext = 0):
    """
    Display an image.
    
    image : numpy array
        2D or a 3D image array
        
    Returns
    -------
    None
    """
    import matplotlib.pylab as plt
    
    plt.figure(figsize = (10,7))

    if image.ndim == 2:
        plt.imshow(image, origin = 'lower', vmin = image.min(), vmax = image.max())
        plt.colorbar()
        plt.axis("off")
        plt.show()
    elif image.ndim == 3:
        nframes = image.shape[0]
        if ext < nframes:
            plt.imshow(image[ext,:,:], origin = 'lower', vmin = image[ext,:,:].min(), vmax = image[ext,:,:].max())
            plt.colorbar()
            plt.axis("off")
            plt.show()
        else:
            print "Warning: Image extension is out of bound. Stopping."
            
            
def imhist(image, ext = 0, bins = 50):
    """
    Plot histogram of image pixel values.
    
    Parameters
    ----------
    image : numpy array
        2D or 3D image array
        
    ext : int
        FITS image extension
        
    bins : int
        Number of histogram bins
        
    Returns
    -------
    None
    """
    import matplotlib.pylab as plt
    
    plt.figure(figsize = (8,8))
    plt.ylabel("Frequency")
    plt.xlabel("Pixel value (ADU)")
        
    if len(image.shape) == 2:
        ny, nx = image.shape
        plt.hist(image.reshape(ny*nx), bins = bins)
        plt.show()
    elif len(image.shape) == 3:
        nframes, ny, nx = image.shape
        if ext < nframes:
            plt.hist(image[ext,:,:].reshape(ny*nx), bins = bins)
            plt.show()
        else:
            print "Warning: Image extension is out of bound. Stopping."
            
    return
    
    
def imcombine(img, combine = "average", nframes = 100):
    """
    Combine image frames of CHIMERA instrument image cube.
    
    Parameters
    ----------
    img : numpy array
        CHIMERA instrument image cube
        
    combine : string
        Combine type - average, sum or median
    
    nframes : int
        Number of frames to combine
        
    Returns
    -------
    comb_image : numpy array
        Combined image numpy array
    """
    ny, nx = img.shape[1], img.shape[2]
    
    if combine == "sum":
        comb_data = np.int16(np.sum(img[:nframes], axis = 0, dtype = np.float32))
    elif combine == "median":
        comb_data = np.zeros([ny, nx], dtype = np.float32)
        np.median(img[:nframes], axis = 0, out = comb_data)
    else:
        comb_data = np.mean(img[:nframes], axis = 0, dtype = np.float32)
        
    return comb_data