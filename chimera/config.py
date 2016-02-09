
import os, json
import chimera


class Config:
    def __init__(self, fname = "chimera.conf"):
        self.conf_fname = os.path.join(chimera.pkg_dir, fname)
        
    def dump(self):
        """
        Generate CHIMERA configuration json file.
        
        Parameters
        ----------
        
        Returns
        -------
        
        """
        cfg = {
            "Detector" : 
            { 
                "8680": 
                {
                    "CONVENTIONAL": 
                    {
                                        "1.0": 
                                        {
                                                "1x": (3.87, 8.51),
                                                "2.4x": (1.6, 6.74),
                                                "4.9x": (0.72, 6.23)
                                        },
                                                    
                                        "3.0": 
                                        {
                                                "1x": (10.12, 14.07),
                                                "2.4x": (4.2, 11.17),
                                                "4.9x": (1.89, 10)
                                        }
                                                        
                    },
                                                  
                                                  
                    "EMGAIN": 
                    {
                                        "1.0" : 
                                        {
                                                "1x": (18.97, 32.25),
                                                "2.4x": (7.61, 19.41),
                                                "4.9x": (3.47, 16.31)
                                        },
                                                
                                        "3.0" : 
                                        {
                                                "1x": (46.56, 54.01),
                                                "2.4x": (19.82, 33.3),
                                                "4.9x": (8.84, 26.25)
                                        },
                                                
                                        "5.0" : 
                                        {
                                                "1x": (46.49, 70.66),
                                                "2.4x": (19.53, 45.11),
                                                "4.9x": (8.9, 35.87)
                                        },
                                                
                                        "10.0" : 
                                        {
                                                "2.4x": (22.45, 52.98),
                                                "4.9x": (10.43, 45.37),
                                        }                                                                                                                                                                                                                                                                                                                                                
                    } 
                                        }, 
                                
                "8325": 
                {
                    "CONVENTIONAL": 
                    {
                                        "1.0": 
                                        {
                                                "1x": (3.98, 8.64),
                                                "2.5x": (1.6, 6.75),
                                                "5.1x": (0.72, 6.23)
                                        },
                                        
                                        "3.0": 
                                        {
                                                "1x": (10.45, 14.42),
                                                "2.5x": (4.14, 10.97),
                                                "5.1x": (1.89, 10.24)
                                        }
                    },
                    
                    "EMGAIN": 
                    {
                                        "1.0" : 
                                        {
                                                "1": (19.73, 34.13),
                                                "2.5x": (7.88, 20.49),
                                                "5.1x": (3.54, 16.99)
                                        },                                      
                                                
                                        "3.0" : 
                                        {
                                                "1": (48.23, 54.5),
                                                "2.5x": (19.77, 33.41),
                                                "5.1x": (9.04, 27.84)
                                        },
                                                
                                        "5.0" : 
                                        {
                                                "1": (50.66, 77.0),
                                                "2.5x": (20.46, 48.08),
                                                "5.1x": (8.7, 35.5)
                                        },
                                                
                                       "10.0" : 
                                       {
                                                "2.5x": (22.44, 53.63),
                                                "5.1x": (11.3, 52.55),
                                       }                                       
                    }
                }
            },
            
            "Phot" :
            {
               "fwhmpsf": 6.0,
               "sigma": 10.0,
               "exposure": "EXPTIME",
               "calgorithm": "centroid",
               "cbox" : 8,
               "maxshift": 5,  
               "salgorithm": "median",
               "annulus": 14,
               "dannulus": 16,
               "apertures": 12,
               "zmag": 27.11
            }
        }
                        
         
        # Dump the configuration to json output file
        with open(self.conf_fname, "w") as fd:
            json.dump(cfg, fd)        
                              
        return
        
        
    def load(self):
        """
        Load the json configuration file.
        
        Parameters
        ----------
        
        Returns
        -------
        cfg : dict
            CHIMERA configuration dictionary
        """
        with open(self.conf_fname, "r") as fd:
            config = json.load(fd)
            
        return config