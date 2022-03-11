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
from pIMOS.read import lisst as read_lisst

font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

class LISST(xrwrap):
    
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
        self.update_attribute('title', 'Measured data from a LISST Data Logger')
        self.update_attribute('source', 'Seabird Data Logger') # Could be more specific.
        self.update_attribute('history', '')
        self.update_attribute('comment', '')

    def read(self, method):
        
        print('Reading info')
        fields, units, dtype, bins_lower, bins_median = read_lisst.load_LISST_info(self.folder)
        print('Reading csv')
        df = pd.read_csv(self.fullpath, header=None)
        print('Converting time')
        time_series = read_lisst.convert_LISST_time(df)
        print('Adding flags')
        df = read_lisst.add_LISST_flags(df)
        print('Parsing csv')
        ds = read_lisst.parse_LISST_csv(df, fields, units, dtype, bins_lower)

        self.ds = ds

        # This is using ALL Will's code - need to discuss c.f. compliance with Will

        # Associate QC Flags - will nead Need to discuss this with Will. 
        print('Done')
        
        return self.ds
