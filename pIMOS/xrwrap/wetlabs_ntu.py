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
import datetime
import os 
import pdb

import importlib 

# import zutils.xrwrap as xrwrap
import pIMOS.xrwrap.xrwrap as xrwrap
import pIMOS.xrwrap.pimoswrap as pimoswrap
import pIMOS.utils.file as zfile


class_attrs = {
            'title': 'Measured data from a WetLABS data logger',
            'source': 'pIMOS', # Could be more specific.
            'process_level': 1
        }


def from_log(filename):
    """
    Spin up from serial log file
    """
        
    skiprows=0
    line = 'n'
    with open(filename) as fop:
        while not 'records to read' in line:
            skiprows+=1
            line = fop.readline()
        pass

    df = pd.read_csv(filename, 
                     skiprows=skiprows, 
                     delim_whitespace=True, 
                     names=['date', 'time', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9'],
                     encoding='latin-1')

    df.time = df.date + ' ' + df.time
    df = df.drop(['date'], axis=1)
    df.time = pd.to_datetime(df.time,dayfirst = True, format = "%m/%d/%y %H:%M:%S")

    df.drop(df.tail(1).index,inplace=True) 

    ds = df.set_index('time').to_xarray()

    ds.attrs['raw_file_name']      = os.path.split(filename)[1]
    ds.attrs['raw_file_directory'] = os.path.split(filename)[0]

    rr = WETLABS_NTU(ds)

    return rr, ds


def from_netcdf(infile):
    """
    Pass straight to the main xrwrap from_netcdf method.
    """
   
    classhandler = WETLABS_NTU

    rr, ds = xrwrap._from_netcdf(infile, classhandler)
    
    return rr, ds

##########################
# Actual xarray wrap #####
##########################
class WETLABS_NTU(pimoswrap.pimoswrap):

    class_attrs = class_attrs
    
    def __init__(self, ds, verbose=True):
        self.verbose = verbose

        if self.verbose:
            print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)
        self.update_attributes_with_dict(class_attrs)
        self.enforce_these_attrs(class_attrs)
        
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

            if 'turbidity' in calibration_data.keys():
                SF = calibration_data['turbidity']['SF'] 
                DC = calibration_data['turbidity']['DC']
                column = int(calibration_data['turbidity']['dev_file_column'] - 2)
                self._calibrate_turbidity(SF, DC, column)

            if 'chlorophyll' in calibration_data.keys():
                SF = calibration_data['chlorophyll']['SF'] 
                DC = calibration_data['chlorophyll']['DC']
                column = int(calibration_data['chlorophyll']['dev_file_column'] - 2)
                self._calibrate_chlorophyll(SF, DC, column)
        
            print('Device file calibration complete.')

    def _calibrate_turbidity(self, SF, DC, column):
        """
        Calibrate turbidity with scale factor and dark counts
        """

        self._calibrate_SFDC(SF, DC, column, 'turbidity', 'seawater_turbidity', 'NTU')
        self.associate_qc_flag('turbidity', 'turbidity')

    def _calibrate_chlorophyll(self, SF, DC, column):
        """
        Calibrate chlorophyll with scale factor and dark counts
        """

        self._calibrate_SFDC(SF, DC, column, 'chlorophyll', 'seawater_chlorophyll_concentration', 'mg/L')
        self.associate_qc_flag('chlorophyll', 'chlorophyll')

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

    
    def flag_maxvalues(self, varname, devfile, tolerance=0.015, max_counts=4130, verbose=None):
        if verbose is None:
            verbose = self.verbose

        unit, model, serial, date_created, calibration_data = parse_device_file(devfile)

        SF = calibration_data[varname]['SF']
        var_max = max_counts * SF

        # Check for max values
        max_flag = self.ds[varname] >= var_max * (1 - tolerance)
        print(np.sum(max_flag.values))

        # Update QC flag
        self.update_qc_flag_logical(self.ds[varname].attrs['qc_variable'], 'time', max_flag, 1, verbose=False)

        
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
            calibration_data['chlorophyll'] = {'SF': float(words[1]), 'DC':int(words[2]), 'dev_file_column':int(i+1)}

        elif param == 'ntu':
            calibration_data['turbidity'] = {'SF': float(words[1]), 'DC':int(words[2]), 'dev_file_column': int(i+1)}
            
        else:
            raise(DeviceFileFormatError('Unrecognised parameter "{}"'.format(param)))

    return unit, model, serial, date_created, calibration_data

class DeviceFileFormatError(Exception):
    pass

class DeviceFileReadError(Exception):
    pass

#%%
