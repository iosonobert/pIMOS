from setuptools import setup, find_packages
from setuptools.dist import Distribution
import os

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

# Get the version info We do this to avoid importing __init__, which
# depends on other packages that may not yet be installed.
base_dir = os.path.abspath(os.path.dirname(__file__))
version = {}
with open(base_dir + "/pIMOS/_version.py") as fp:
    exec(fp.read(), version)

setup(name='pIMOS',
      version=version['__version__'],
      url='https://github.com/iosonobert/pIMOS/',
      description='field data processing',
      author='Andrew Zulberti',
      author_email='andrew.zulberti@gmail.com',
      packages=['pIMOS.utils', 
                'pIMOS.xrwrap',
                'pIMOS.read',],
      # packages=find_packages(),
      install_requires=['numpy==1.24.2',
                        'matplotlib>=3.6.3', 
                        'xarray==2023.1.0',
                        'dolfyn>=1.2.0',
                        'turbo_lance'],
      license='unlicensed to all but author',
      include_package_data=True,
      distclass=BinaryDistribution,
    )
