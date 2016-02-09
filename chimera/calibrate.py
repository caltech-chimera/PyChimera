
from __future__ import division

import numpy as np


def masterbias(image, skip = 0):
    """
    Create an average master bias frame.
    
    Parameters
    ----------
    image : numpy array
        3D bias array from CHIMERA
        
    skip : int
        Number of frames of skip from start. Default is 0 (no frame is skipped).
        
    Returns
    -------
    avg_bias : numpy array
        2D average bias image
    """
    if image.ndim == 3:
        avg_bias = np.mean(image[skip:,:,:], axis = 0, dtype = np.float32)
        return avg_bias
    else:
        print "MASTERBIAS: Only 3D image arrays are supported"
        
    return
    
    	
def masterflat(flat_image, bias_image, threshold = 0.8):	
    """
    Generate normalized field.
    
    Parameters
    ----------
    flat_image : numpy array
        3D flat field image
        
    bias_image : numpy array
        2D master bias image
        
    threshold : float
        Normalized threshold used to normalize the flat field
        
    Returns
    -------
    norm_flat : numpy array
        Normalized flat field image
    """
    # Check if flat fields are 3D array and bias frames is 2D
    if flat_image.ndim == 3 and bias_image.ndim == 2:
        # Subtract bias dc from flat fields
        flat_image -= bias_image
    
        # Normalize each frame by mean value
        nframes, ny, nx = flat_image.shape
        norm_flat_data = np.zeros([nframes, ny, nx])
        for i in range(nframes):
            norm_flat_data[i,:,:] = flat_image[i,:,:] / np.mean(flat_image[i,:,:])
    
        # Median normalized flat field
        median_flat = np.median(norm_flat_data, axis = 0)
    
        mask = median_flat > threshold
        median_flat /= np.mean(median_flat[mask])
        return median_flat
    else:
        print "MASTERFLAT: Only 3D flat field and 2D master bias supported."
        return
    


def debias(sci_image, bias_image):
    """
    De-bias science frames.
    
    Parameters
    ----------
    sci_image : numpy array
        Raw science image
        
    bias_image : numpy array
        Master bias image
    
    Returns
    -------
    debias_sci_image : numpy array
        Debias science image
    """
    if bias_image.ndim == 2 and (sci_image.ndim == 2 or sci_image.ndim == 3):
        debias_sci_image = sci_image - bias_image
        return debias_sci_image
    else:
        print "DEBIAS: Only 2D master bias and 2D/3D science images are supported."
        return
    
    
def imreduce(sci_image, bias_image, flat_image):
    """
    De-bias and flat field the science frames as well as generate an average
    science frame (if the tracking is poor, streaked star image with be 
    generated.)
    
    Parameters
    ----------
    sci_image : numpy array
        2D or 3D science image array
        
    bias_image : numpy array
        2D master bias image array
        
    flat_image : numpy array
        2D normalized flat field image array
        
    Returns
    -------
    sci_red_image : numpy array
        2D or 3D reduced image array
        
    avg_sci_image : numpy array
        2D average science image array
        
    Note: Average science image is only returned if the input science image 
    array is 3D otherwise only reduced science image is returned.        
    """
    # De-bias science frames
    debias_sci_image = sci_image - bias_image
    
    # Flat-field the science frames
    sci_red_image = debias_sci_image/flat_image
    sci_red_image = np.asarray(sci_red_image, dtype = np.float32)
    
    # Average the frames
    if sci_image.ndim == 3:
        avg_sci_image = np.mean(sci_image, axis = 0, dtype = np.float32)
        return sci_red_image, avg_sci_image
    else:
        return sci_red_image