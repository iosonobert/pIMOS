#!/usr/bin/env python

from setuptools import setup, find_packages
from setuptools.dist import Distribution

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

		
setup(name='turbo_tools',
      version='1.0.0',
      description='field turbulence calcs',
      author='Andrew Zulberti',
      author_email='andrew.zulberti@gmail.com',
      #packages=['zutils'],
      packages=find_packages(),
      install_requires=['numpy','matplotlib', 'xarray'],
      license='unlicensed to all but author',
      include_package_data=True,
      distclass=BinaryDistribution,
    )
