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

# import zutils.xrwrap as xrwrap
import pIMOS.xrwrap.xrwrap as xrwrap
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

    rr.add_variable_attributes()

    return rr, ds

def from_cnvfile(filename):
        
    ds = read_sbdctd.read(filename)

    ds.attrs['raw_file_name']      = os.path.split(filename)[1]
    ds.attrs['raw_file_directory'] = os.path.split(filename)[0]

    rr = SOLANDER_CTD(ds)
        
    # This is toward process level 1 stuff
    rr.add_variable_attributes()
    
    return rr, ds
        
##########################
# Actual xarray wrap #####
##########################
class SOLANDER_CTD(xrwrap.xrwrap):
    
    class_attrs = {
            'title': 'Measured data from a profiling Seabird CTD',
            'source': 'pIMOS'
        }

    def __init__(self, ds):
        
        print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)
        
        self.enforce_these_attrs(self.class_attrs)

        pass

    
    def add_variable_attributes(self):
        """
        Adds attributes to variables and associates QC flags to variables. 

        Not sure if this is exhaustive. May need continual improvement. 
        """

        if 'depSM' in self.ds:
            self.ds = self.ds.rename(name_dict={'depSM': 'Depth'})

        if 'prDM' in self.ds:
            self.ds = self.ds.rename(name_dict={'prDM': 'Pressure'})
        
        if 't090C' in self.ds:
            self.ds = self.ds.rename(name_dict={'t090C': 'Temperature'})
            
        if 't190C' in self.ds:
            self.ds = self.ds.rename(name_dict={'t190C': 'Temperature2'})
        
        if 'c0Sperm' in self.ds:
            self.ds = self.ds.rename(name_dict={'c0Sperm': 'Conductivity'})

        if 'c1Sperm' in self.ds:
            self.ds = self.ds.rename(name_dict={'c1Sperm': 'Conductivity2'})

        if 'density00' in self.ds:
            self.ds = self.ds.rename(name_dict={'density00': 'Density'})
            
        if 'flECO_AFL' in self.ds:
            self.ds = self.ds.rename(name_dict={'flECO_AFL': 'Fluorescence'})
            
        if 'par' in self.ds:
            self.ds = self.ds.rename(name_dict={'par': 'PAR'})
            
        if 'sal00' in self.ds:
            self.ds = self.ds.rename(name_dict={'sal00': 'Salinity'})
        
        if 'sal11' in self.ds:
            self.ds = self.ds.rename(name_dict={'sal11': 'Salinity2'})
        
        if 'CStarTr0' in self.ds:
            self.ds = self.ds.rename(name_dict={'CStarTr0': 'Beam_transmission'})
        
        if 'CStarAt0' in self.ds:
            self.ds = self.ds.rename(name_dict={'CStarAt0': 'Beam_attenuation'})
        
        if 'altM' in self.ds:
            self.ds = self.ds.rename(name_dict={'altM': 'Altimeter'})
        
        if 'accM' in self.ds:
            self.ds = self.ds.rename(name_dict={'accM': 'Acceleration'})
        
        if 'timeS' in self.ds:
            self.ds = self.ds.assign_coords({'time': self.ds.timeS})
            self.ds = self.ds.drop('timeS')
            
        ds = self.ds

        for v in ['Temperature', 'Temperature2']:
            if v in ds.data_vars:
                ds[v].attrs['long_name'] = 'seawater_temperature'
                ds[v].attrs['standard_name'] = 'seawater_temperature'
                ds[v].attrs['units'] = 'deg'
                ds[v].attrs['cf_compliant'] = 1
                self.associate_qc_flag(v, v)

        if 'Pressure' in ds.data_vars:
            ds['Pressure'].attrs['standard_name'] = 'sea_water_pressure'
            ds['Pressure'].attrs['long_name'] = '"Sea water pressure" is the pressure that exists in the medium of sea water. It includes the pressure due to overlying sea water, sea ice, air and any other medium that may be present. For sea water pressure excluding the pressure due to overlying media other than sea water, the standard name sea_water_pressure_due_to_sea_water should be used.'
            ds['Pressure'].attrs['units'] = 'dbar'
            ds['Pressure'].attrs['cf_compliant'] = 1
            self.associate_qc_flag('Pressure', 'Pressure')

        for v in ['Conductivity', 'Conductivity2']:
            if v in ds.data_vars:
                ds[v].attrs['long_name'] = 'sea_water_electrical_conductivity'
                ds[v].attrs['standard_name'] = 'sea_water_electrical_conductivity'
                ds[v].attrs['units'] = 'S m-1'
                ds[v].attrs['cf_compliant'] = 1
                self.associate_qc_flag(v, v)
        
        if 'Density' in ds.data_vars:
            ds['Conductivity'].attrs['long_name'] = 'sea_water_density'
            ds['Conductivity'].attrs['standard_name'] = 'sea_water_density'
            ds['Conductivity'].attrs['units'] = 'kg m-3'
            ds['Conductivity'].attrs['cf_compliant'] = 1
            self.associate_qc_flag('Density', 'Density')

        for v in ['Salinity', 'Salinity2']:
            if v in ds.data_vars:
                ds[v].attrs['long_name'] = 'sea_water_salinity'
                ds[v].attrs['standard_name'] = 'sea_water_salinity'
                ds[v].attrs['units'] = 'PSU'
                ds[v].attrs['cf_compliant'] = 1
                self.associate_qc_flag(v, v)

        if 'Depth' in ds.data_vars:
            ds['Depth'].attrs['long_name'] = 'pressure_sensor_depth_below_sea_surface'
            ds['Depth'].attrs['standard_name'] = 'pressure_sensor_depth_below_sea_surface'
            ds['Depth'].attrs['units'] = 'm'
            ds['Depth'].attrs['cf_compliant'] = 1
            self.associate_qc_flag('Depth', 'Pressure')
            
        for var in ds:
            if not 'cf_compliant' in ds[var].attrs:
                print('Setting {} to be non cf_compliant'.format(var))
                ds[var].attrs['cf_compliant'] = 0
                
