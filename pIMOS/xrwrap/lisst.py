"""
Created on Mon Jul 30 15:17:42 2018

@author: 20270917
"""
#%%

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
from matplotlib.dates import num2date, date2num
import matplotlib
import scipy.signal as signal
import pdb
import datetime
import os

# import zutils.xrwrap as xrwrap
import pIMOS.xrwrap.xrwrap as xrwrap

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
    
    class_attrs = {
            'title': 'Measured data from a LISST Data Logger',
            'source': 'pIMOS' # Could be more specific.
        }

    def __init__(self, ds):
        
        print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        
        self.enforce_these_attrs(self.class_attrs)
