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
from pIMOS.read import SEABIRD_37_39_56 as read_sbd

font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)


def from_asc(filename):
    """
    Spin up from *.asc file
    """
    ds = read_sbd.parse_seabird_asc(filename)

    ds.attrs['raw_file_name']      = os.path.split(filename)[1]
    ds.attrs['raw_file_directory'] = os.path.split(filename)[0]
        
    rr = SEABIRD_37_39_56(ds)
        
    # This is toward process level 1 stuff
    rr.add_variable_attributes()
    
    return rr, ds

def from_cnv(filename):
    """
    Spin up from *.asc file
    """
        
    ds = read_sbd.parse_seabird_cnv(filename)

    ds.attrs['raw_file_name']      = os.path.split(filename)[1]
    ds.attrs['raw_file_directory'] = os.path.split(filename)[0]
        
    rr = SEABIRD_37_39_56(ds)
    
    # This is toward process level 1 stuff
    rr.add_variable_attributes()

    return rr, ds

def from_netcdf(infile):
    """
    Pass straight to the main xrwrap from_netcdf method.
    """
   
    classhandler = SEABIRD_37_39_56

    rr, ds = xrwrap._from_netcdf(infile, classhandler)
    
    return rr, ds

##########################
# Actual xarray wrap #####
##########################
class SEABIRD_37_39_56(xrwrap.xrwrap):
    
    # folder = ''
    # file = ''

    so = None
    eo = None
    
    job = ''
    trip = ''
   
    def __init__(self, ds):
        
        print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        class_attrs = {
            'title': 'Measured data from a Seabird Data Logger',
            'source': 'Seabird Data Logger' # Could be more specific.
        }
        self.enforce_these_attrs(class_attrs)


    def add_variable_attributes(self):
        """
        Adds attributes to variables and associates QC flags to variables. 
        """
        ds = self.ds

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
