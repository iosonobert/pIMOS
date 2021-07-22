# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:17:42 2018

@author: 20270917
"""
#%%
import matplotlib.pyplot as plt
import pandas
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib
from windrose import WindroseAxes
import datetime
import os 
import pdb

import importlib 

from zutils.xrwrap import xrwrap
import zutils.stats as zstats
import zutils.file as zfile

class WETLABS_NTU(xrwrap):

    def __init__(self, infile, driver='pandas', attributes={}):
        
        self.parse_infile(infile)
        
        self.driver_name = driver

        spam_loader = importlib.find_loader(driver)
        if spam_loader is None:
            raise(Exception('No module found for driver {}.'.format(driver)))

        self.driver = importlib.import_module(driver)

        if driver.lower() in ['pandas']: 
            self.read_raw_pd(attributes=attributes)
        elif driver.lower() == 'xarray':
            self.load(self.folder, self.file_)
        else:
            raise(Exception('{} is not a valid driver'.format(driver)))

        self.update_global_attrs()
    
    def update_global_attrs(self):
        """
        Each wrapper should overload this function
        """

        # CF Compliance
        print('Updating attributes function of the class.')
        self.update_attribute('title', 'Measured data from a TDRI ADCP read from .PD0 files')
        self.update_attribute('institution', 'UWA')
        self.update_attribute('source', 'TDRI ADCP [Workhorse, Quartermaster, or Longranger]')
        self.update_attribute('history', '')
        self.update_attribute('references', '')
        self.update_attribute('comment', '')
        
    def read_raw_pd(self, attributes):
        
        df = pd.read_csv(self.fullpath, skiprows=1, delim_whitespace=True, names=['date', 'time', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'])

        df.time = df.date + ' ' + df.time
        df = df.drop(['date'], axis=1)
        df.time = pd.to_datetime(df.time,dayfirst = True, format = "%m/%d/%y %H:%M:%S")

        df.drop(df.tail(1).index,inplace=True) 

        ds = df.set_index('time').to_xarray()

        self.ds = ds
        
    def _calibrate_device(self, device_file):
        """
        Calibrate using device file. Read the SF, DC and columns.
        """

        unit, model, serial, date_created, calibration_data = parse_device_file(device_file)

        print('Wetlab model:' + model)
        print('Wetlab ' + model + 'serial number:' + serial)
        print('Date device file created:' + date_created)
        
        if len(calibration_data) == 0:
            print('File contains no calibration information')
        else:
            print('File contains  calibration information for:')
            for key in calibration_data.keys():
                print('     ' + key)

            if 'Turbidity' in calibration_data.keys():
                SF = calibration_data['Turbidity']['SF'] 
                DC = calibration_data['Turbidity']['DC']
                column = int(calibration_data['Turbidity']['dev_file_column'] - 2)
                self._calibrate_turbidity(SF, DC, column)

            if 'Chlorophyll' in calibration_data.keys():
                SF = calibration_data['Chlorophyll']['SF'] 
                DC = calibration_data['Chlorophyll']['DC']
                column = int(calibration_data['Chlorophyll']['dev_file_column'] - 2)
                self._calibrate_chlorophyll(SF, DC, column)
        
            print('Device file calibration complete.')

    def _calibrate_turbidity(self, SF, DC, column):
        """
        Calibrate turbidity with scale factor and dark counts
        """

        self._calibrate_SFDC(SF, DC, column, 'turbidity', 'seawater_turbidity', 'NTU')

    def _calibrate_chlorophyll(self, SF, DC, column):
        """
        Calibrate chlorophyll with scale factor and dark counts
        """

        self._calibrate_SFDC(SF, DC, column, 'chlorophyll', 'seawater_chlorophyll_concentration', 'mg/L')

    def _calibrate_SFDC(self, SF, DC, column, param_name, param_long_name, param_units):
        """
        Calibrate parameter with scale factor and dark counts
        """

        ds = self.ds
        
        if not type(column) in [int, np.int32]:
            raise(Exception("Column number must be an integer"))
            
        column_name = 'c{:0}'.format(column)
        if not column_name in ds.data_vars.keys():
            raise(Exception("Invalid column number"))
            
        print(column_name)
        
        ds[param_name] = xr.DataArray(SF*(ds[column_name].values-DC), dims={'time':ds.time})
        ds[param_name].attrs['units'] = param_units
        ds[param_name].attrs['long_name'] = param_long_name
        ds[param_name].attrs['standard_name'] = param_long_name
        
    def _calc_burst_mean(self, parameter='turbidity', cutoff=None):
        
        """
        cutoff: the time in seconds that distinguishes bursts
        """
        
        ds = self.ds
        
        if cutoff is None: 
            cutoff = np.mean(np.diff(ds['time']))/10
        else:
            cutoff = np.timedelta64(cutoff, 's')
        print('Splitting bursts using a cutoff of {}'.format(cutoff))

        dt = np.diff(ds['time'])
        burst_ind = np.where(np.diff(ds['time'].values)>cutoff)[0]+1

        burst_ind = np.concatenate((np.array([0]), burst_ind, np.array([len(ds['time'].values)])))

        s = burst_ind.shape

        burst_mean = np.zeros((s[0]-1,))
        burst_median = burst_mean.copy()
        burst_median_mad = burst_mean.copy()
        burst_std = burst_mean.copy()
        burst_mean_time = burst_mean.copy().tolist()

        for i in np.arange(0, len(burst_mean)):
            ind = np.arange(burst_ind[i], burst_ind[i+1])

            burst_mean_time[i]   = np.mean(ds['time'][ind]).values
            burst_mean[i]   = np.mean(ds[parameter][ind])
            burst_median[i] = np.median(ds[parameter][ind])
            burst_std[i] = np.var(ds[parameter][ind])**(1/2)
            burst_median_mad[i] = np.max(np.abs(ds[parameter][ind]-burst_median[i])) # Max Absolute Deviation

            pass

        burst_mean_time = np.array(burst_mean_time)

        name='burst_mean_'+parameter
        ds[name] = xr.DataArray( burst_mean,\
            dims=('burst_time'),\
            name=name,\
            attrs = ds[parameter].attrs,\
            coords = {'burst_time': burst_mean_time}
        )

        name='burst_median_'+parameter
        ds[name] = xr.DataArray( burst_median,\
            dims=('burst_time'),\
            name=name,\
            attrs = ds[parameter].attrs,\
            coords = {'burst_time': burst_mean_time}
        )

        name='burst_std_'+parameter
        ds[name] = xr.DataArray( burst_std,\
            dims=('burst_time'),\
            name=name,\
            attrs = ds[parameter].attrs,\
            coords = {'burst_time': burst_mean_time}
        )
        
        name='burst_median_mad_'+parameter
        ds[name] = xr.DataArray( burst_median_mad,\
            dims=('burst_time'),\
            name=name,\
            attrs = ds[parameter].attrs,\
            coords = {'burst_time': burst_mean_time}
        )

def parse_device_file(dev_file, verbose=False):

    outputs = _parse_device_file(dev_file, verbose=False)

    return outputs

def _parse_device_file(dev_file, verbose=False):
    calibration_data = {}
    with open(dev_file) as f:

        lines = f.readlines()
        if verbose:
            for line in lines:
                print(line)

    # Line 0 expecting something like "ECO FLNTUSB-1835"
    ECO, unit = lines[0].split()
    model, serial = unit.split('-')

    # Line 1 expecting something like "Created on: 	mm/dd/yy"
    Created_on, date_created = lines[1].split(':')
    if not Created_on.lower()=='created on':
        raise(DeviceFileFormatError)
    date_created = date_created.strip()

    # Line 2 expecting "\n"

    # Line 3 expecting "COLUMNS=7n\n"
    COLUMNS, columns = lines[3].split('=')
    if not COLUMNS.lower()=='columns':
        raise(DeviceFileFormatError)
    columns = int(columns)

    for i in np.arange(0, columns):
        words = lines[4+i].split()

        param, d_file_param_col = words[0].lower().split('=')

        if param == 'n/u':
            continue

        elif param == 'chl':
            calibration_data['Chlorophyll'] = {'SF': float(words[1]), 'DC':int(words[2]), 'dev_file_column':int(i+1)}

        elif param == 'ntu':
            calibration_data['Turbidity'] = {'SF': float(words[1]), 'DC':int(words[2]), 'dev_file_column': int(i+1)}
            
        else:
            raise(DeviceFileFormatError('Unrecognised parameter "{}"'.format(param)))

    return unit, model, serial, date_created, calibration_data

class DeviceFileFormatError(Exception):
    pass

class DeviceFileReadError(Exception):
    pass

#%%
