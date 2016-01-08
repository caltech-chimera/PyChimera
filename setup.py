import ez_setup
ez_setup.use_setuptools()

from setuptools import find_packages, setup

setup(
      name="CHIMERA",
      version = "0.1",
      packages=find_packages(),
      
      # PROPER uses numpy and scipy
      install_requires = ['numpy>=1.8', 'scipy>=0.14', 'pyfits>=3.0', 'pyraf>=2.1', 'photutils>=0.2'],
      
      package_data = {
        # If any package contains *.txt, *.rst or *.fits files, include them:
        '': ['*.txt', '*.rst', '*.fits']
      },
      
      # Metadata for upload to PyPI
      author="Navtej Singh Saini",
      author_email = "navtej@astro.caltech.edu",
      description="CHIMERA instrument data processing pipeline",
      license = "BSD",
      platforms=["OSX", "Linux", "Unix"],
      url="https://github.com/navtej/CHIMERA",
) 
