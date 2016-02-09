import os.path as _osp

pkg_dir = _osp.abspath(_osp.dirname(__file__))


from .animate import animate
from .imutil import imfill, imshow, imhist, imcombine
from .fitsutil import fitsread, fitswrite, fitshead, fitscombine
from .calibrate import masterbias, masterflat, imreduce
from .aperphot import Aperphot
from .centroid import recenter

# Generate CHIMERA configuration file
from chimera import config
cfg = config.Config()
cfg.dump()