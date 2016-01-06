
from __future__ import division

import chimera
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def animate(image, nframes, vmin, vmax, scale, cmap, fps, outtype):
    """ 
    Generate animation from CHIMERA image cube.
    
    Parameters
    ----------
    image : string
        CHIMERA FITS image name
        
    nframes : int
        Number of frames to include in animation
        
    vmin, vmax : float
        Minimum and maximum pixel value
        
    scale : string
        Scale the images? (linear or log)
        
    cmap : string
        Matlab colormap
    
    fps : int
        Frames per second for mp4 animation
        
    outtype : string
        Save animation as gif or mp4
    
    Returns
    -------
    anim_fname : string
        Name of the animation file name.
    """
    
    print "ANIMATE: Animate CHIMERA 3D image cubes"
    
    nframes = int(nframes)
    vmin = int(vmin)
    vmax = int(vmax)
    fps = int(fps)
    
    # Create a matplotlib figure
    fig = plt.figure(figsize = (12,10))
    plt.axis("off")
    
    # Read FITS image
    img_data = chimera.fitsread(image)        
    
    # Display each frame to create an animation
    if nframes > img_data.shape[0]:
        nframes = img_data.shape[0]
    
    print "  Generating animation"
         
    ims = []
    for i in range(0, nframes):
        if scale == 'linear':
            im = plt.imshow(img_data[i,:,:], origin = "lower", cmap = plt.get_cmap(cmap), vmin = vmin, vmax = vmax)
        else:
            img_log = np.log10(img_data[i,:,:])
            im = plt.imshow(img_log, origin = "lower", cmap = plt.get_cmap(cmap), vmin = np.min(img_log), vmax = np.max(img_log))
        ims.append([im])

    # Animate the frames
    ani = animation.ArtistAnimation(fig, ims, interval = 50, blit = True, repeat_delay = 10)

    # Save the animation
    if outtype == "gif":
        anim_fname = image.replace(".fits", ".gif")
        ani.save(anim_fname, dpi = 100, writer = "imagemagick")
    else:
        anim_fname = image.replace(".fits", ".mp4")
        ani.save(anim_fname, fps = fps, dpi = 100, extra_args=['-vcodec', 'libx264'])
                
    return anim_fname