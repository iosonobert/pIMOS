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

from zutils.xr import xrwrap
import pIMOS.read.SEABIRD_CTD as read_sbdctd

font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

class SOLANDER_CTD(xrwrap):
    folder = ''
    file = ''
    so = None
    eo = None
    
    job = ''
    trip = ''
    
    def __init__(self, folder, file_, attributes={}, model=None, method='csv'):
        
        self.folder = folder
        self.file_ = file_
        self.attributes = attributes

        if method.lower() in ['csv']: 
            self.read(method=method)
        elif method.lower() == 'nc':
            raise(Exception('Must implement NC reading'))
        elif method.lower() == 'xr':
            raise(Exception('Must implement XR wrapping'))
        else:
            raise(Exception('Unknown method'))

        self.update_global_attrs()
        self.update_attributes_with_dict(attributes)

    def update_global_attrs(self):
        """
        Each wrapper should overload this function
        """

        print('Setting default attributes of the class.')
        self.update_attribute('title', 'Measured data from a profiling Seabird CTD')
        self.update_attribute('institution', 'UWA')
        self.update_attribute('source', 'Profiling CTD from A.I.M.S. R.V. Solander') # Could be more specific.
        self.update_attribute('history', '')
        self.update_attribute('references', 'CF Conventions version 1.7')
        self.update_attribute('comment', '')

    def read(self, method):
        
        self.ds = read_sbdctd.read(self.fullpath)

        print('Done')
        
        return self.ds
