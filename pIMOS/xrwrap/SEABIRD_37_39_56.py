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
from pIMOS.read import SEABIRD_37_39_56 as read_sbd

font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

class SEABIRD_37_39_56(xrwrap):
    
    folder = ''
    file = ''
    so = None
    eo = None
    
    job = ''
    trip = ''
    
    def __init__(self, folder, file_, attributes={}, model=None, method=None):
        
        self.folder = folder
        self.file_ = file_
        self.attributes = attributes

        # Make sure a model is given
        if model is None:
            raise(Exception('Must specify a model.'))
        if not model.lower() in ['sbe56', 'sbe39', 'sbe37']:
            raise(Exception('{} is not a recognised logger.'.format(model)))

        # If no method given, run the default method for each model. The default was decided by AZ based in Apr 2021 by testing each method.
        # There is no known reason why all methods should not work. 
        if method is None:
            if model.lower() in ['sbe39']:
                method = 'asc'
            if model.lower() in ['sbe37', 'sbe56']:
                method = 'cnv'

        if method.lower() in ['asc', 'cnv']: 
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
        self.update_attribute('title', 'Measured data from a Seabird Data Logger')
        # self.update_attribute('institution', 'O2 Metocean')
        self.update_attribute('source', 'Seabird Data Logger') # Could be more specific.
        self.update_attribute('history', '')
        # self.update_attribute('references', 'O2 Metocean QAQC Conventions; CF Conventions version 1.7')
        self.update_attribute('comment', '')

    def read_data(self, method):
        
        if method=='asc':
            ds = read_sbd.parse_seabird_asc(self.fullpath)
        if method=='cnv':
            ds = read_sbd.parse_seabird_cnv(self.fullpath)

        return ds

    def read(self, method):
        
        ds = self.read_data(method)
        self.ds = ds

        # Ultimately I'll need to modify this to be more dynamic. This is why we call it crude_readers.
        if 'Temperature' in ds.data_vars:
            ds['Temperature'].attrs['long_name'] = 'seawater_temperature'
            ds['Temperature'].attrs['standard_name'] = 'seawater_temperature'
            ds['Temperature'].attrs['units'] = 'deg'
            self.associate_qc_flag('Temperature', 'Temperature')

        if 'Pressure' in ds.data_vars:
            ds['Pressure'].attrs['standard_name'] = 'sea_water_pressure'
            ds['Pressure'].attrs['long_name'] = '"Sea water pressure" is the pressure that exists in the medium of sea water. It includes the pressure due to overlying sea water, sea ice, air and any other medium that may be present. For sea water pressure excluding the pressure due to overlying media other than sea water, the standard name sea_water_pressure_due_to_sea_water should be used.'
            ds['Pressure'].attrs['units'] = 'dbar'
            self.associate_qc_flag('Pressure', 'Pressure')

        if 'Conductivity' in ds.data_vars:
            ds['Conductivity'].attrs['long_name'] = 'sea_water_electrical_conductivity'
            ds['Conductivity'].attrs['standard_name'] = 'sea_water_electrical_conductivity'
            ds['Conductivity'].attrs['units'] = 'S m-1'
            self.associate_qc_flag('Conductivity', 'Conductivity')

        # ds['Sea pressure'].attrs['standard_name'] = 'sea_water_pressure_due_to_sea_water'
        # ds['Sea pressure'].attrs['long_name'] = 'The pressure that exists in the medium of sea water due to overlying sea water. Excludes the pressure due to sea ice, air and any other medium that may be present. For sea water pressure including the pressure due to overlying media other than sea water, the standard name sea_water_pressure should be used.'
        # ds['Sea pressure'].attrs['units'] = 'dbar'
        # self.associate_qc_flag('Sea pressure', 'Pressure')

        if 'Depth' in ds.data_vars:
            ds['Depth'].attrs['long_name'] = 'pressure_sensor_depth_below_sea_surface'
            ds['Depth'].attrs['standard_name'] = 'pressure_sensor_depth_below_sea_surface'
            ds['Depth'].attrs['units'] = 'm'
            self.associate_qc_flag('Depth', 'Pressure')

        return self.ds

    def export(self, final=False):

        if final:
            self.ds.attrs['Disclaimer'] = self.disclaimer
            outname = outname + 'finalised'
            if not final_folder is None:
                folder = final_folder

        self.ds.close() # Force close
        self.ds.attrs = self._attrs 

        self.ds.to_dataframe().to_csv('{folder}//{file_}.csv'.format(folder=self.folder, file_=self.file_))
        self.ds.to_netcdf(path='{folder}//{file_}.nc'.format(folder=self.folder, file_=self.file_))

        self.ds.close() # Force close

        try: # Again this double DS nonsense is a bit of a nightmare. Need a better solution to this. 
            self.ds_wave.close() # Force close
            self.ds.attrs = self._attrs 
            self.ds_wave.to_dataframe().to_csv('{folder}//{file_}_wave.csv'.format(folder=self.folder, file_=self.file_))
            self.ds_wave.to_netcdf(path='{folder}//{file_}_wave.nc'.format(folder=self.folder, file_=self.file_))
            self.ds_wave.close() # Force close
        except:
            pass


# %% Wave code below from Mike Cuttler

def disperk(f, h):

    # % DISPERK   Linear dispersion relation
    # %
    # % usage:   k  = disper(w,h,g)     or
    # %
    # %          k  = wave number             (2 * pi / wave length)
    # %          w  = wave angular frequency  (2 * pi / wave period)
    # %          h  = water depth
    # %          g  = gravitational acceleration constant, optional (DEFAULT 9.81)
    # %
    # %          absolute error in k*h < 5.0e-16 for all k*h
    # %

    # %          programmer: G. Klopman, Delft Hydraulics, 6 Dec 1994

    w = 2*np.pi*f

    g = 9.81
    eps = 2.22e-16
    
    w2 = (w**2)* h/ g
    q  = w2/ (1 - np.exp (-(w2**(5/4))))**(2/5)

    for j in np.arange(0, 3):
        thq     = np.tanh(q)
        thq2    = 1 - thq**2
        a       = (1 - q * thq) * thq2
        b       = thq + q * thq2
        c       = q * thq - w2
        arg     = (b**2) - 4 * a * c
        arg     = (-b + np.sqrt(arg)) / (2 * a+eps)
        iq      = np.where(abs(a*c) < 1.0e-8 * (b**2))

        arg[iq] = - c[iq] / b[iq]
        q       = q + arg

    k = np.sign(w)*q/h
    
    return k

def get_spec(x, correct=False, site_depth=None, hab=0.5, f_co=0.4):
    """
    f_co = cutoff frequency. Just set everything above to zero
    """
    f, Pxx_raw = signal.csd(x, x, fs=2.0, nperseg=256, noverlap=128)
    k = disperk(f, site_depth)
    
    if correct:
        correction = np.cosh(k*site_depth)/np.cosh(k*hab)
        correction = correction**2
    
        correction[np.where(f>f_co)] = 0

        Pxx = Pxx_raw*correction
        
        # Pxx[np.where(f>f_co)] = 0
        
    return f, k, Pxx, Pxx_raw, correction
