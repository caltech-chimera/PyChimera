
import pyfits
import chimera
import numpy as np
from pyraf import iraf
from datetime import datetime, timedelta

    
MONTHS = {"Jan": 1, "Feb": 2, "March": 3, "April": 4, "May": 5, "June": 6, "July": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12}


class Aperphot:
    def __init__(self, sci_file, coords):
        self.sci_file = sci_file
        self.coords = coords
        
        # Set header keyword parameters
        self.setkeywords()
        
        # Set IRAF parameters and load packages
        self.setiraf()
        self.setparams()


    def setiraf(self):
        """
        Set IRAF global parameters and load DAOPHOT package for aperture
        photometry.
        
        Parameters
        ----------
        
        Returns
        -------
        None
        """
        iraf.prcacheOff()
        iraf.set(writepars=0)
    
        # Load IRAF packages
        iraf.noao(_doprint = 0)
        iraf.noao.digiphot(_doprint = 0)
        iraf.noao.digiphot.daophot(_doprint = 0)
    
        return
        
                        
    def setparams(self):
        """
        Set datapars, centerpars, fitskypars and photpars.
        
        Parameteres
        -----------
        
        Returns
        -------
        None
        """  
        # Set parameters for daophot 
        self.fwhmpsf = 6.0
        self.sigma = 10.0
        self.readnoise = 10
        self.epadu = 1.89
        self.exposure = 'EXPTIME'

        self.calgorithm = 'centroid'
        self.cbox = 8
        self.maxshift = 5

        self.salgorithm = 'median' 
        self.annulus = 14 
        self.dannulus = 16

        self.apertures = 12.0
        self.zmag = 27.11
        
        # Set parameters for phot
        self.method = "exact"
        self.inner_radius = 14
        self.outer_radius = 30

        return
        

    def setkeywords(self):
        """
        Set FITS image header keyword parameters.
        
        Parameters
        ----------
        
        Returns
        -------
        None
        """
        header = pyfits.getheader(self.sci_file, ignore_missing_end = True)
        self.nx = header["NAXIS1"]
        self.ny = header["NAXIS2"]
        self.nframes = header["NAXIS3"]
        self.kintime = header["KINCYCTI"]
        utcstart = header["UTCSTART"]
        self.utcstart = self.parser(utcstart)
        
                        
    def parser(self, utcstart):
        """
        Datetime parser for CHIMERA UTCSTART header keyword.
        
        Parameters
        ----------
        utcstart : string
            Datetime for start of frame (in UTC)
            
        Returns
        -------
        dt : datetime struct
            Datetime structure
        """
        month, date, year, time = utcstart.split("-")
        month = MONTHS[month]
        date = int(date)
        year = int(year)
        
        hour, minu, sec = time.split(":")
        hour = int(hour)
        minu = int(minu)
        sec, ms = sec.split(".")
        sec = int(sec)
        ms = int(ms) * 1000
        
        dt = datetime(year, month, date, hour, minu, sec, ms)
        
        return dt
    

    def addtime(self, secs):
        """
        Add time in seconds to UTC datetime.
        
        Parameters
        ----------
        secs : float
            Time to add to UTC in seconds.
            
        Returns
        -------
        dt : datetime structure
        """    
        td = timedelta(0, secs)

        return self.utcstart + td
        
            
    def daocog(self, tolerance = 0.01):
        """
        Curve of growth to determine nominal aperture for photometry using DAOPHOT.
        
        Parameters
        ----------
        tolerance : float
            Magnitude difference tolerance between different apertures
        
        Returns
        -------
        aperture : float
            Nominal aperture radius for photmetry
        """
        # Randomly peform curve of growth on 5 frames
        framenum = np.random.randint(1, self.nframes, 5)
        
        apertures = np.linspace(2,20,19)
        
        # Iterate through the frames and determine nominal aperture
        nom_aper = np.zeros(5, dtype = np.float32)
        cnt = 0
        for val in framenum:
            outfile = self.sci_file.replace(".fits", "." + str(val) + ".cog.phot.1")
            iraf.delete(outfile)
            self.daophot(val, self.coords, outfile, apertures = "2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20")
            mags = iraf.pdump(outfile, "mag", "yes", Stdout = 1)
            mags_arr = np.array(mags[1].split(),dtype = np.float32)
            mags_diff = np.diff(mags_arr)
            idx = np.where(np.abs(mags_diff) < 0.01)
            if len(idx[0]) != 0:
                nom_aper[cnt] = apertures[idx[0][0]]
            else:
                nom_aper[cnt] = 12.0
            cnt += 1
            iraf.delete(outfile)
            
        return np.median(nom_aper)
        

    def cog(self, window_size, method, tolerance = 0.01):
        """
        Curve of growth to determine nominal aperture for photometry using 
        astropy photutils.
        
        Parameters
        ----------
        tolerance : float
            Magnitude difference tolerance between different apertures
        
        Returns
        -------
        aperture : float
            Nominal aperture radius for photmetry
        """
        # Aperture values in pixels
        apertures = np.array([2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
        naper = apertures.shape[0]
        
        # Randomly peform curve of growth on 5 frames
        framenum = np.random.randint(1, self.nframes, 5)
        
        apertures = np.linspace(2,20,19)
        
        # Read input image and star position
        image = chimera.fitsread(self.sci_file)
        pos = np.loadtxt(self.coords, ndmin = 2)
        
        # Iterate through the frames and determine nominal aperture
        nom_aper = np.zeros(5, dtype = np.float32)
        cnt = 0
        for val in framenum:
            mags_arr = np.zeros(len(apertures))
            objpos = chimera.recenter(image[val,:,:], pos, window_size, method)
            for i in range(naper):
                flux = self.phot(image[val,:,:], objpos, aper = apertures[i])
                mags_arr[i] = -2.5 * np.log10(flux['flux'])
            mags_diff = np.diff(mags_arr)
            idx = np.where(np.abs(mags_diff) < 0.01)
            if len(idx[0]) != 0:
                nom_aper[cnt] = apertures[idx[0][0]]
            else:
                nom_aper[cnt] = 12.0
            cnt += 1
            
        return np.median(nom_aper)
        
                
    def daophot(self, framenum, coords, outfile, apertures, verbose = "no"):
        """
        Aperture photometry of stars in the coords file using IRAF PHOT routine.
        
        Parameters
        ----------
        framenum : int
            Frame number in the image cube to perform aperture photometry on.
            
        coords : string
            Text file with coordinate of the stars
            
        outfile : string
            Text file to which photometry results are written
        
        Returns
        -------
        None
        """
        iraf.delete(outfile)
        iraf.phot(image = self.sci_file + "[,," + str(framenum) + "]", coords = coords, output = outfile, fwhmpsf = self.fwhmpsf, sigma = self.sigma, readnoise = self.readnoise, epadu = self.epadu, exposure = self.exposure, calgorithm = self.calgorithm, cbox = self.cbox, maxshift = self.maxshift, salgorithm = self.salgorithm, annulus = self.annulus, dannulus = self.dannulus, apertures = apertures, zmag = self.zmag, interactive = "no", verify = "no", verbose = verbose)        
        
        return
        
        
    def phot(self, image, objpos, aper):
        """
        Aperture photometry using Astropy's photutils.
        
        Parameters
        ----------
        image : numpy array
            2D image array
            
        objpos : list of tuple
            Object poistions as list of tuples
            
        aper : float
            Aperture radius in pixels
         
        Returns 
        -------
        phot_table : astropy table
             Output table with stellar photometry   
        """
        try:
            from astropy.table import hstack
            from photutils import aperture_photometry, CircularAnnulus, CircularAperture
        except ImportError:
            pass
    
        apertures = CircularAperture(objpos, r = aper) 
        annulus_apertures = CircularAnnulus(objpos, r_in = self.inner_radius, r_out = self.outer_radius)
        
        rawflux_table = aperture_photometry(image, apertures = apertures, method = self.method)
        bkgflux_table = aperture_photometry(image, apertures = annulus_apertures, method = self.method)
        phot_table = hstack([rawflux_table, bkgflux_table], table_names = ["raw", "bkg"])
        
        bkg = phot_table["aperture_sum_bkg"] / annulus_apertures.area()
        phot_table["msky"] = bkg
        phot_table["area"] = apertures.area()
        phot_table["nsky"] = annulus_apertures.area()
                
        bkg_sum = bkg * apertures.area()
        final_sum = phot_table["aperture_sum_raw"] - bkg_sum
        phot_table["flux"] = final_sum
        
        return phot_table