#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Script to load raw LISST deployment data from CSV files to NetCDF files
 - Set for randomly shaped model files

l_lisst_to_nc.py

### CODE AUTHOR ###
 * William Edge

### CODE UPDATE ###
 * 20181011: Adaption of code from profiling script
 * 20190530: Update to handle fieldtrip input 
'''


### IMPORT PACKAGES ###
import os
import glob
import pandas as pd
import numpy as np
import sys

# fieldtrip = 'RS2019'
# sn = '2031' # 2009
# bg = 'D3'

fieldtrip = 'KISSME2017'
sn = '2031' # 2009
bg = 'bottle24'
mooring = 'SP250-L'

datatype = 'Lander'

# Particle shape model
model = 'rs'

# Specify work locations
cur_d = os.path.abspath(os.curdir)
p_loc = cur_d.find('PhD')
i_loc = len(cur_d) - p_loc

workDir = os.path.abspath(os.curdir)[:-(i_loc-3)]

sys.path.insert(0, os.path.join(workDir, 'pl'))
from fielddata import extra
from fielddata import lisst as ls



### CODE INPUT ###
##def main():
if 1==1:
    
    # Specify work locations

    # Define work directories and find LISST data
    rawLoc = os.path.join(workDir,'data','field',fieldtrip,'raw_data',mooring,'LISST')       

    # Source and save folders
    rawDir = os.path.join(rawLoc, sn, bg)
    saveDir = rawDir.replace('raw_data','processed_data')
    # Check is save dir exists and make if not
    if not os.path.isdir(saveDir):
        os.makedirs(saveDir)


    # Find file for particle shape model
    if model=='rs':
        # Find all Randomly shaped CSV files
        rsFileSearch = os.path.join(rawDir, '*rs.csv')
        lis_files = glob.glob(rsFileSearch)
    else:
        # Find all spherical model files
        rsFileSearch = os.path.join(rawDir, '*.csv')
        lis_files = glob.glob(rsFileSearch)
        rem_rs = [str.find('rs.csv') for str in lis_files]
        for ind, p in enumerate(rem_rs):
            if p != -1:
                lis_files.pop(ind)
        
    # Load lisst device info
    fields, units, dtype, bins_lower, bins_median = ls.load_LISST_info(workDir)

##    # Read Fieldtrip CTD logbook
##    logFileName = fieldtrip + '_CTD_Logbook.xlsx' 
##    logbook = os.path.join(workDir,'fieldwork', fieldtrip, logFileName)
##    df_check = pd.read_excel(logbook)


    ## Loop through file list
    for files in lis_files:

        # Read file
        df = pd.read_csv(files, header=None)

        # Convert LISST time to days since 1970
        time_series = ls.convert_LISST_time(df)

        # Set nc filename
        nc_file, nc_full = extra.set_NC_filename(files, saveDir)        

        # Add LISST flags as variables
        new_col = len(df.columns)
        df = ls.add_LISST_flags(df)

        # Add bad flag when depth =-10
        # df.loc[df[40] < 0, [new_col]] = 1

        # Add bad flag when depth > 300 
##        depth_diff = df[40].diff()
        # df.loc[df[40] > 300, [new_col]] = 1

        # Create xarray dataset from subset dataframe
        ds = ls.parse_LISST_csv(df, fields, units, dtype, bins_lower)

        # Add deployment attributes
##        ds.attrs['FieldTripID'] = 'KISSME2017'
##        ds.attrs['Instrument_SerialNumber'] = '2031'
##        ds.attrs['SiteID'] = 'SP250-L'
##        ds.attrs['Filename'] = files
##        ds.attrs['TimeSwitchOn'] = extra.dayssince_2_datestr(time_series[0])
##        ds.attrs['TimeFirstWet'] = '2017-04-01 02:49:12'
##        ds.attrs['TimeFirstInPos'] = '2017-04-01 03:37:36'
##        ds.attrs['TimeFirstGoodData'] = '2017-04-01 03:39:36'
##        ds.attrs['TimeLastGoodData'] = '2017-05-22 03:11:38'
##        ds.attrs['TimeLastInPos'] = '2017-05-22 03:12:00'
##        ds.attrs['TimeOnDeck'] = '2017-05-22 03:33:36'
##        ds.attrs['TimeSwitchOff'] = extra.dayssince_2_datestr(time_series[-1])
##        ds.attrs['InstrumentDepth'] = 272

        ds.attrs['FieldTripID'] = fieldtrip
        ds.attrs['Instrument_SerialNumber'] = sn
        ds.attrs['SiteID'] = 'L150'
        ds.attrs['Filename'] = files
        ds.attrs['TimeSwitchOn'] = extra.dayssince_2_datestr(time_series[0])
        ds.attrs['TimeFirstWet'] = '2019-03-05 00:59:44'
        ds.attrs['TimeFirstInPos'] = '2019-03-05 01:20:24'
        ds.attrs['TimeFirstGoodData'] = '2019-03-05 01:22:24'
        ds.attrs['TimeLastGoodData'] = '2019-04-24 23:12:04'
        ds.attrs['TimeLastInPos'] = '2019-04-24 23:12:04'
        ds.attrs['TimeOnDeck'] = '2019-04-24 23:22:24'
        ds.attrs['TimeSwitchOff'] = extra.dayssince_2_datestr(time_series[-1])
        ds.attrs['InstrumentDepth'] = 147

        ds.attrs['DataStatus'] = 'Auto and manual QAQC flags added'
        ds.attrs['ModifiedBy'] = 'William Edge'        

        # Create NetCDF
        ls.create_LISST_netCDF(nc_file, saveDir, ds, datatype, sn)

##    return None
	

##if __name__ == '__main__':
##	main()


##nc_fid = Dataset(nc_full, 'r', format='NETCDF4')
##ncdump(nc_fid, verb=True)
##nc_fid.close()
