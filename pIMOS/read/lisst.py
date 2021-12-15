#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Package of  LISST related data processing modules

lisst.py

### CODE AUTHOR ###
 * William Edge

### CODE UPDATE ###
 * 20181009: Creation of code
'''


import os
import glob
import pandas as pd
import datetime as dT
import numpy as np
from netCDF4 import Dataset
import time
import xarray as xr
import sys
from pIMOS.utils.othertime import DaysSince


### Objects ###

## Load LISST Output Info
def load_LISST_device(file):

    """ Function to load LISST-200X device info 
        and output lists of field names, units, and data specifiers
        :param: file: LISST output info csv
        :return: fields: list of output fields from instrument by column
        :return: units: list of units for above fields
        :return: dtype: list of data types for above fields (for SQL) 
    """
            
    # Read file
    df = pd.read_csv(file, encoding = "utf-8")

    # Drop separated time info (y,m,d,H,M,S)
    df = df.drop(df.index[43:49], axis=0)

    # Convert dataframe to 3 lists
    fields = df.loc[:,'lisst_csv_fields'].tolist()
    units = df.loc[:,'lisst_units'].tolist()
    dtype = df.loc[:,'data_specifier'].tolist()

    return fields, units, dtype


## Load LISST size bins
def load_LISST_bins(file):

    """ Function to load LISST-200X sediment bin sizes
        :param: file: LISST sediment size bins csv
        :return: bins_lower: list of lower bound of each bin in micrometres
    """
    
    # Read file
    df = pd.read_csv(file, encoding = "utf-8")

    # Convert lowwer size limit bins to list
    bins_lower = df.loc[:,'sizebins_lower'].tolist()
    bins_median = df.loc[:,'sizebins_median'].tolist()

    return bins_lower, bins_median


## Load all LISST info in one go
def load_LISST_info(workDir):

    """ Call both LISST loading modules
        :param: None
        :return: all variables from load functions
    """

    # Load LISST device file data (field names, units, netCDF data specifiers)
    lisst_device_file = os.path.join(workDir, 'LISST_OutputInfo.csv')
    lisst_sizebins = os.path.join(workDir, 'LISST_Sediment_SizeBins.csv')
    fields, units, dtype = load_LISST_device(lisst_device_file)

    # Add info for QAQC flag columns
    fields.append('BadFlag')
    fields.append('WarningFlag')
    units.append('NA')
    units.append('NA')
    dtype.append('i1')
    dtype.append('i1')

    bins_lower, bins_median = load_LISST_bins(lisst_sizebins)

    return fields, units, dtype, bins_lower, bins_median


## Change LISST time to Days Since
def convert_LISST_time(df):

    """ Function to convert time from LISST dataframe to  
        Days Since 1970 using mraysons othertime module
        :param: dataframe of LISST file
        :return: timeout: timeseries in days since 1970        
    """

    #### Must be a better way to convert time ####
    # Pull out time columns 
    yr = df.loc[:,42].tolist()
    mon = df.loc[:,43].tolist()
    dy = df.loc[:,44].tolist()
    hr = df.loc[:,45].tolist()
    mn = df.loc[:,46].tolist()
    sc = df.loc[:,47].tolist()
    micro = np.zeros((len(df),), dtype=int).tolist()

    # Loop to append each datetime to list
    pytime = []
    for a, b, c, d, e, f, g in zip(yr, mon, dy ,hr, mn, sc, micro):
            pt = dT.datetime(int(a),int(b),int(c),int(d),\
                             int(e),int(f),int(g))    
            pytime.append(pt)
            
    # Convert time to days since 1970
    timeout = DaysSince(pytime, basetime = dT.datetime(1970,1,1))  

    return timeout


## Add bad and warning flag columns to dataframe
def add_LISST_flags(df):

    """ Function add flags as per LISST user manual       
    """

    # Create flag columns
    new_col = len(df.columns)
    bad_flag = pd.DataFrame({new_col: np.zeros(len(df))})
    warning_flag = pd.DataFrame({(new_col+1): np.zeros(len(df))})

    # Add new columns to existing df
    df = pd.concat([df, bad_flag], axis=1)
    df = pd.concat([df, warning_flag], axis=1)
    
    # Add conditions to bad flag column
    df.loc[df[59] > 0.995, [new_col]] = 1
    df.loc[df[59] < 0.1, [new_col]] = 1
    
    # Add conditions to warning flag column
    df.loc[df[59] > 0.98, [(new_col+1)]] = 1
    df.loc[df[59] < 0.3, [(new_col+1)]] = 1

    # Remove bad flag from warning column
    df.loc[df[new_col] == 1, [(new_col+1)]] = 0

    return df


## Find where file should be split
def calc_LISST_splits(checkDF):

    """ Function to check if LISST file spans multiple CTD profiles
        and find time to split if required
        :param: checkDF: data from CTD logbook excel file
        :return: splitStart: (list of) start of new file
        :return: splitEnd: (list of) end of new file
    """

    # Get split times for dataset
    splitStart = []
    splitEnd = []
    num_splits = list(range(len(checkDF)))

    # Loop thorugh number of required splits
    for count in num_splits:

        # Calculate split time unless it is last in batch (which goes to end)
        if count != max(num_splits):
                
            # Find average between timeout and next time in
            hr_out = checkDF.iloc[count,3]
            day_out = checkDF.iloc[count,0]
            hr_in = checkDF.iloc[count+1,2]
            day_in = checkDF.iloc[count+1,0]

            # Adjust if profile crosses date line
            prev_in = checkDF.iloc[count,3]
            if int(hr_out[0:2]) < int(prev_in[0:2]):
                print(hr_out)
                print(hr_in)
                dT_out = dT.datetime.strptime((str(day_out)[0:10] + hr_out), '%Y-%m-%d%H:%M')\
                                 + dT.timedelta(days=1)
                print(dT_out)
            else:
                dT_out = dT.datetime.strptime((str(day_out)[0:10] + hr_out), '%Y-%m-%d%H:%M')
                    
            dT_in = dT.datetime.strptime((str(day_in)[0:10] + hr_in), '%Y-%m-%d%H:%M')

            # Find middle of earlier profile end and later profile start
            sp_mean = dT_out + (dT_in - dT_out)/2
            print('File requires split at %s' % str(sp_mean))
                 
        # Compile splits 
        if count==0:
            # Split from start of file
            splitStart.append(0)
            splitEnd.append(DaysSince(sp_mean, basetime = dT.datetime(1970,1,1)))
        elif count==max(num_splits):
            # Split to end of file
            splitStart.append(splitEnd[count-1])
            splitEnd.append(999999)
        else:
            # Split between two identified breaks in profiling
            splitStart.append(splitEnd[count-1])
            splitEnd.append(DaysSince(sp_mean, basetime = dT.datetime(1970,1,1)))

    return splitStart, splitEnd


## Read and convert LISST dataframe to a dataset 
def parse_LISST_csv(df, nfields, nunits, ndtypes, bins_lower):
	
    """ Function to read df and convert to an XArray dataset with
        variable name, unit, & dataype
        :param: df: LISST raw dataframe from CSV
        :param: LISST device info variables
        :return: ds: dataset of LISST data
    """

    # Convert time
    timeout = convert_LISST_time(df)

    # Drop old time columns from dataframe
    df_drop = df.drop([42,43,44,45,46,47], axis=1) 
    df_drop.columns = range(df_drop.shape[1])

    # Create conc dataframe 
    df_conc = df.loc[:,0:35]

    # Create an output dataset
    ds = xr.Dataset()

    # Loop through df columns corresponding to device file field names
    ind = -1	
    for field, unit, dtype in zip(nfields, nunits, ndtypes):
            
        # Set normal attributes
        normattrs = {}
        normattrs['dims'] = ('time',)
        normattrs['name'] = field
        normattrs['attrs'] = {'units':unit, 'name':field}
        normattrs['coords'] = {'time':timeout}

        # Set conc attributes
        concattrs = {}
        concattrs['dims'] = ('time', 'bins_lower',)
        concattrs['name'] = 'conc'
        concattrs['attrs'] = {'units':unit, 'name':'conc', 'data_specifier':dtype}
        concattrs['coords'] = {'time':timeout, 'bins_lower':bins_lower}

        # First variable is time
        if ind==-1:
                data = timeout
                allattrs = normattrs
                
        # Second variable is conc dataframe
        elif ind==0:
                data = df_conc
                ind = ind + 35
                allattrs = concattrs
                
        # All other variables come from dataframe in order  	
        elif ind>35:
                data = df_drop.loc[:,ind].tolist()
                allattrs = normattrs

        V = xr.DataArray( \
                data, \
                dims = allattrs['dims'], \
                name = allattrs['name'],\
                attrs = allattrs['attrs'],\
                coords = allattrs['coords']
                )     

        ds.update({str(field):V})
        ind+=1

    return ds


## Create new name for split file
def name_split_file(nc_file, count, model):

    """ Function to rename new split file with alpha appendage
        :param: nc_file: original NC filename
        :param: count: position of split file (int)
        :param: model: data procesing model (spherical or rs)
        :return: file_new: new NC filename
    """
    
    # Set new split filename starting with A, B, etc
    if model=='rs':
            file_new = nc_file[:-3] + '-' + chr(count + ord('A')) + '_rs'
    else:
            file_new = nc_file + '-' + chr(count + ord('A'))
            
    print('Raw file %s being split into %s' % (nc_file, file_new))

    return file_new


## Start new file
def initiate_new_nc(nc_full, nc_file):

    """ Function to initiate new netCDF file and delete  
        old file with same name if exists
        :param: nc_full: full path to new NC file inc name
        :param: nc_file: new NC file name only
        :return: rootgrp: new NC file handle
        :return: rootcheck: 
    """

    # Delete existing netCDF file (if exist)
    if os.path.exists(nc_full):
        os.remove(nc_full)
        print('Existing NetCDF file deleted')

    # Create file
    rootgrp = Dataset(nc_full, "w", format="NETCDF4")
    rootgrp.name = nc_file
    rootcheck = True
    print('Initialising file %s' % nc_file)

    return rootgrp, rootcheck


## Write new NC file
def create_LISST_netCDF(nc_file, save_path, data, datatype, sn):

    """ Write processed NC file and save
        :param: nc_file: new NC file name only
        :param: save_path: full path to save NC file
        :param: data: dataset variable
        :param: datatype: string indicating data collection type 
        :param: sn: instrument serial number
    """

    # Set error check variable to false
    rootcheck = False

    # Try create file and delete if errors
    try:
        # Set approriate filename extension
        if nc_file[-3] != '.nc':
                nc_file += '.nc'

        # Set nc filename
        nc_full = os.path.join(save_path, nc_file)

        # Create file template
        rootgrp, rootcheck = initiate_new_nc(nc_full, nc_file)

        # Add netCDF global attributes
        rootgrp.description = datatype + ' LISST Data from SN: ' + sn         
        rootgrp.created = 'Created ' + time.ctime(time.time())
        rootgrp.created_by = 'William Edge'
##        rootgrp.QAQC_level = 'None'
        rootgrp.close()

        # Add dimensions to time variable
        data.time.attrs['units'] = "days since 1970-01-01 00:00:00.0"

        # Write dataset to new netcdf file
        data.to_netcdf(nc_full, 'a', 'NETCDF4')

        print("File %s saved successfully\n" % nc_file)

    # Remove NC file if exception occurs
    except Exception as e: 
        print(e)
        print(nc_file)
        if rootcheck:
            rootgrp.close()
            os.remove(nc_full)
            print("""Error saving file %s 
                            File deleted""" % nc_file)

    return None


