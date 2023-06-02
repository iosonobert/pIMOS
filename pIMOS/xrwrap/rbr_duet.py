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
# from windrose import WindroseAxes
import scipy.signal as signal
import pdb
import datetime
import os

# import zutils.xrwrap as xrwrap
from pyrsktools import RSK

import pIMOS.xrwrap.xrwrap as xrwrap

font = {'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)


def from_rsk(filename):
    """
    Spin up from *.rsk file
    """
        
    
    # Instantiate an RSK class object, passing the path to an RSK file
    rsk = RSK(filename)
    # Open the RSK file. Metadata is read here
    print('Opening file....')
    rsk.open()
    # Read, process, view, or export data here
    print('Reading file....')
    rsk.readdata()

    longNames = [c.longName for c in rsk.channels]
    longNames

    ti = [i for i in [0, 1] if longNames[i] == 'temperature'][0] + 1
    pi = [i for i in [0, 1] if longNames[i] == 'pressure'][0] + 1

    print('Really slow loops here.....')
    time = [d[0] for d in rsk.data]
    temperature = [d[ti] for d in rsk.data]
    pressure = [d[pi] for d in rsk.data]
    print('.....done')

    # Close the RSK file
    rsk.close()
    print('File closed.')


    ds = xr.Dataset({       'Pressure': ('time', pressure),
                            'Temperature': ('time', temperature)},
                            coords={'time': time})
        
    rr = RBR_DUET(ds)
    
    # This is toward process level 1 stuff
    rr.add_variable_attributes()

    return rr, ds

def from_netcdf(infile):
    """
    Pass straight to the main xrwrap from_netcdf method.
    """
   
    classhandler = RBR_DUET

    rr, ds = xrwrap._from_netcdf(infile, classhandler)
    
    return rr, ds

##########################
# Actual xarray wrap #####
##########################
class RBR_DUET(xrwrap.xrwrap):

    class_attrs = {
            'title': 'Measured data from an RBR Duet Data Logger',
            'source': 'RBR Duet Data Logger' # Could be more specific.
        }

    def __init__(self, ds):
        
        print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        self.enforce_these_attrs(self.class_attrs)


    def add_variable_attributes(self):
        """
        Adds attributes to variables and associates QC flags to variables. 
        """
        ds = self.ds

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

        if 'Depth' in ds.data_vars:
            ds['Depth'].attrs['long_name'] = 'pressure_sensor_depth_below_sea_surface'
            ds['Depth'].attrs['standard_name'] = 'pressure_sensor_depth_below_sea_surface'
            ds['Depth'].attrs['units'] = 'm'
            self.associate_qc_flag('Depth', 'Pressure')


        return self.ds
