"""
Created on Mon Jul 30 15:17:42 2018

@author: 20270917
"""
#%%


#%%

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
from matplotlib.dates import num2date, date2num
import matplotlib
from windrose import WindroseAxes
import scipy.signal as signal
import pdb
import datetime
import os

import zutils.xrwrap as xrwrap
import pIMOS.read.SEABIRD_CTD as read_sbdctd

font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

def from_netcdf(infile):
    """
    Pass straight to the main xrwrap load method.
    """

    folder, file = xrwrap.parse_infile(infile)

    ds = xr.open_dataset(os.path.join(folder, file))

    ds.attrs['last_load_file_name']      = file
    ds.attrs['last_load_directory']      = folder

    print(ds)
    
    rr = SOLANDER_CTD(ds)

    return rr, ds

def from_cnvfile(filename):
        
    ds = read_sbdctd.read(filename)

    ds.attrs['raw_file_name']      = os.path.split(filename)[1]
    ds.attrs['raw_file_directory'] = os.path.split(filename)[0]

    rr = SOLANDER_CTD(ds)
        
    return rr, ds
        
##########################
# Actual xarray wrap #####
##########################
class SOLANDER_CTD(xrwrap.xrwrap):
    folder = ''
    file = ''
    so = None
    eo = None
    
    job = ''
    trip = ''
    
    # def __init__(self, folder, file_, attributes={}, model=None, method='csv'):
        
    #     self.folder = folder
    #     self.file_ = file_
    #     self.attributes = attributes

    #     if method.lower() in ['csv']: 
    #         self.read(method=method)
    #     elif method.lower() == 'nc':
    #         raise(Exception('Must implement NC reading'))
    #     elif method.lower() == 'xr':
    #         raise(Exception('Must implement XR wrapping'))
    #     else:
    #         raise(Exception('Unknown method'))

    #     self.update_global_attrs()
    #     self.update_attributes_with_dict(attributes)

    def __init__(self, ds):
        
        print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        class_attrs = {
            'title': 'Measured data from a profiling Seabird CTD',
            'source': 'Profiling CTD from A.I.M.S. R.V. Solander'
        }
        self.enforce_these_attrs(class_attrs)

        pass

    
