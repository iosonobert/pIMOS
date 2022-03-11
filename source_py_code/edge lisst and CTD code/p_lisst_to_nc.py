#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Script to load raw LISST profiling data from CSV files to NetCDF files
 - Set for randomly shaped model files

p_lisst_to_nc.py

### CODE AUTHOR ###
 * William Edge

### CODE UPDATE ###
 * 20180717: Creation of code
 * 20180821: QAQC from manual added
 * 20190520: Multiple fieldtrip handling added
'''


### IMPORT PACKAGES ###
import os
import glob
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, os.path.abspath(os.curdir)[:-9])
from fielddata import extra
from fielddata import lisst as ls


##fieldtrip = 'RS2019'
##sn = '2031'
##bg = 'D4'
fieldtrip = 'KISSME2017'
sn = '2009'
bg = 'DI2'


### CODE INPUT ###
##def main():
if 1==1:
    
    # Specify work locations
    workDir = os.path.abspath(os.curdir)[:-12]

    # Define work directories and find LISST data
    rawLoc = os.path.join(workDir,'data','field',fieldtrip,'raw_data','profiling',\
                          'calibration_data','LISST')
     
    datatype = 'Profiling'  

    # Source and save folders
    rawDir = os.path.join(rawLoc, sn, bg)
    saveDir = rawDir.replace('raw_data','processed_data')
    # Check is save dir exists and make if not
    if not os.path.isdir(saveDir):
        os.makedirs(saveDir)

    # Particle shape model
    model = 'rs'
    
    if model=='rs':
        # Find all Randomly shaped CSV files
        rsFileSearch = os.path.join(rawDir, '*rs.csv')
        csv_files = glob.glob(rsFileSearch)
        

    fields, units, dtype, bins_lower, bins_median = ls.load_LISST_info(workDir)

    # Read CTD logbook
    logFileName = fieldtrip + '_CTD_Logbook.xlsx' 
    logbook = os.path.join(workDir,'fieldwork',fieldtrip,logFileName)
    df_check = pd.read_excel(logbook)

    ## Loop through file list and append together
    df_master = pd.DataFrame([])
    for files in csv_files:

        # Read file
        df = pd.read_csv(files, header=None)

        # Convert LISST time to days since 1970
        time_series = ls.convert_LISST_time(df)

        # Set nc filename
        nc_file, nc_full = extra.set_NC_filename(files, saveDir)        

        df = ls.add_LISST_flags(df)

        # Check if valid filename in logbook
        if model=='rs':
            df_check_temp = df_check[df_check['LISST Raw Filename'].str.match(nc_file[:-3], na=False)]
        else:
            df_check_temp = df_check[df_check['LISST Raw Filename'].str.match(nc_file, na=False)]      


        ## Check if LISST file needs to be split to match CTD a file
        if len(df_check_temp) > 1:

            # Work out split times
            splitStart, splitEnd = ls.calc_LISST_splits(df_check_temp)

            # Split dataset at splitTimes
            count = 0
            for start, end in zip(splitStart, splitEnd):

                ind = np.where(np.logical_and(time_series>start, time_series<end))
                df_new = df.iloc[ind]
                
                # Create xarray dataset from subset dataframe
                ds_new = ls.parse_LISST_csv(df_new, fields, units, dtype, bins_lower)  
                
                # Name new split
                file_new = ls.name_split_file(nc_file, count, model)

                # Create NetCDF
                ls.create_LISST_netCDF(file_new, saveDir, ds_new, datatype, sn)           
                count += 1

        else:
            # Create NetCDF
            ds_new = ls.parse_LISST_csv(df, fields, units, dtype, bins_lower)  
            ls.create_LISST_netCDF(nc_file, saveDir, ds_new, datatype, sn)           


        df_master = df_master.append(df)
      
    # Create xarray dataset from subset dataframe
    ds = ls.parse_LISST_csv(df_master, fields, units, dtype, bins_lower)  

    # Create NetCDF
    ls.create_LISST_netCDF((fieldtrip + '_allcasts'), saveDir, ds, datatype, sn)


