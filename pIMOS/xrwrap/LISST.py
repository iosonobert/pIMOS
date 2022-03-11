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
from pIMOS.read import lisst as read_lisst

font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

def from_csv(filename):
        
        folder, file =  xrwrap.parse_infile(filename)
        fullpath = os.path.join(folder, file)

        print('Reading info')
        fields, units, dtype, bins_lower, bins_median = read_lisst.load_LISST_info(folder)
        print('Reading csv')
        df = pd.read_csv(fullpath, header=None)
        print('Converting time')
        time_series = read_lisst.convert_LISST_time(df)
        print('Adding flags')
        df = read_lisst.add_LISST_flags(df)
        print('Parsing csv')
        ds = read_lisst.parse_LISST_csv(df, fields, units, dtype, bins_lower)

        ds.attrs['raw_file_name']      = os.path.split(filename)[1]
        ds.attrs['raw_file_directory'] = os.path.split(filename)[0]

        rr = LISST(ds)
        # This is using ALL Will's code - need to discuss c.f. compliance with Will

        # Associate QC Flags - will nead Need to discuss this with Will. 
        print('Done')
        
        return rr, ds

class LISST(xrwrap.xrwrap):
    
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
            'title': 'Measured data from a LISST Data Logger',
            'source': 'LISST Data Logger' # Could be more specific.
        }
        self.enforce_these_attrs(class_attrs)
