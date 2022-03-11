#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Script to convert raw CTD .cnv files to NetCDF4

ctd_2_nc.py

### CODE AUTHOR ###
 * William Edge

### CODE UPDATE ###
 * 20180808: Creation of code
 * 20180822: Split into modules
'''


### IMPORT PACKAGES ###
import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset
import time
import datetime
import sys
sys.path.insert(0, os.path.abspath(os.curdir)[:-9])
from fielddata import extra
from mycurrents import seabird
from soda.utils import othertime


def read_cnv_header(file):

    ''' Function to read text header of a CNV file into
        a single string seperated by newlines
    '''

    header = []
    string = "*END*"

    # Read file
    with open(file, 'rt') as in_file:

        # Append each line until END
        for line in in_file:
            if line.find(string) == -1:
                header.append(line)
            else:
                break

        # Join header together into single string with new lines
        header = "\n".join(header)
    return header


def create_CTD_timeseries(sbe, header):

    ''' Function to create a timeseries for the length of a
        CTD file based on 
    '''

    # Create time variable
    time_start = sbe.datetime
        
    # Find timestep
    t_st = header.find('interval = seconds:')
    time_step = float(header[(t_st+20):(t_st+29)])
    time_step = time_step/24/60/60 # adjust to days not seconds

    # Create time series array from time start
    time_series = np.arange(sbe.nrecords)
    time_series = time_series*time_step

    # Add in basetime from 1970
    time_start_1970 = othertime.DaysSince(time_start)
    time_series = time_series + time_start_1970

    return time_series

        
def get_SBE_units(sbe):

    ''' Function to extract units as string from SBE longnames
    '''

    units = []
    units_convert = {' Voltage 2':'volts',\
                   ' PAR/Irradiance, Biospherical/Licor':'mmol/m^2',\
                   '  0.000e+00':''}

    for longstr in sbe.longnames:
        # Find square brackets
        loc1 = longstr.find('[')
        if loc1>0:
            loc2 = longstr.find(']')
            units.append(longstr[(loc1+1):loc2])
        else:
            # Check convert dictionary
            try:
                units.append(units_convert[longstr])
            except:
                units.append('')

    return units


def create_bottle_array(bpos):

    ''' Function to create array where non-zero entry equals bottle fire
        by number and row, each bottle applied over rows
    '''
    
    pos_diff = np.diff(bpos) # find diff between each element of array
    bot_ind = np.nonzero(pos_diff) # find where diff is non-zero
    print(bot_ind)

    # Create array of zeros
    bpos_new = np.zeros(len(bpos))
    sample_length = 36 # bottle open for 36 scans of CTD (~1.5 seconds at 24 Hz)
##    count = 1 # bottle count for each loop
    

    # Loop through each bottle fire and apply count number to sample rows
    for ind in bot_ind[0]:
            bpos_new[ind:(ind+sample_length)] = bpos[ind+1]
##          count += 1

    return bpos_new


def check_bottle_array(sbe):

    '''Ammend bpos array to indicate where bottle actually fired
    '''

    # Find bottle position array
    if 'bpos' in sbe.names:
        bot_pos = sbe.names.index('bpos')    
        bpos = np.array(sbe.array[:,bot_pos])

        # If a bottle was fired, alter array
        if np.nanmax(bpos) > 0:

            bpos_new = create_bottle_array(bpos)

            sbe.array[:,bot_pos] = bpos_new
    return sbe
    


def sbe_to_dataset(ds, sbe, units):

    '''Convert all arrays in an SBE class to arrays in a dataset, ds,
        with dimensions depSM and units specified by a list of unit strings
    '''

    # Loop through arrays in SBE
    for ind, fields in enumerate(sbe.names):
                
        # Attributes spec with matching list of units
        attrs = {'units':units[ind], 'longname':sbe.longnames[ind].strip()}

        # Specify all variables with the dimension depth
        V = xr.DataArray( \
            sbe.array[:,ind], \
            dims=('depSM',), \
            name=fields,\
            attrs = attrs,\
            coords = {'depSM':sbe.array[:,sbe.names.index('depSM')].tolist()})

        # Update dataset with new array           
        ds.update({str(fields):V})

    return ds


def add_ctd_globals(rootgrp, sbe, ds, nc_file, fieldtrip, maxdepth, sitedepth, numbottles, rawfile):

    '''Add ctd global attributes to new netCDF file
    '''

    rootgrp.CTDFilename = nc_file
    rootgrp.Instruments_Make = 'Seabird'
    rootgrp.Instrument_Model = 'SBE911'
    rootgrp.Instrument_SerialNumber = fieldtrip +'_Profiler'
    rootgrp.FieldTripID = fieldtrip
    rootgrp.AvgLatitude = sbe.lat
    rootgrp.AvgLongitude = sbe.lon
    rootgrp.SiteID = ''
    rootgrp.TimeFirstInPos = sbe.datetime.strftime('%Y-%m-%d %H:%M:%S')
    rootgrp.TimeLastInPos = (sbe.datetime + datetime.timedelta(0,ds.timeS.values.max()))\
                              .strftime('%Y-%m-%d %H:%M:%S')
    rootgrp.TimeZone = 'UTC'
    rootgrp.MaxInstrumentDepth = maxdepth
    rootgrp.SiteDepth = sitedepth
    rootgrp.DepthDatum = 'Sea surface'

    if numbottles != 0:
        rootgrp.BottleSamplesTaken = 'YES'
        rootgrp.NumberBottleSamples = numbottles
        print('Bottle globals added')
    else:
        rootgrp.BottleSamplesTaken = 'NO'
        rootgrp.NumberBottleSamples = '0'
    
    rootgrp.description = 'CTD profiling data from ' + fieldtrip
    rootgrp.TimeModified = time.strftime('%Y-%m-%d %H:%M:%S')
    rootgrp.ModifiedBy = 'William Edge'
    rootgrp.QAQCLevel = 'None'
    rootgrp.close()

    return None
    

def create_ctd_nc(nc_full, saveDir, ds, sbe, rawfile, fieldtrip):
    # rootcheck changes to false after file creation 
    rootcheck = False

    # Try to create new file and delete if error occurs
    try:
        nc_file = extra.get_filename(nc_full)
        # Delete existing netCDF file (if exist)
        if os.path.exists(nc_full):
            os.remove(nc_full)
            print('Existing NetCDF file deleted')
        
        # Create file
        rootgrp = Dataset(nc_full, "w", format="NETCDF4")
        rootcheck = True
        print('Initialising file %s.nc' % nc_file)
        
        # Add netCDF global attributes
        maxdepth = ds.depSM.values.max()
        i_max = ds.depSM.values.argmax()
        sitedepth = ds.depSM.values.max() + ds.altM[i_max]
        numbottles = ds.nbf.values.max()
        add_ctd_globals(rootgrp, sbe, ds, nc_file, fieldtrip, maxdepth, sitedepth, numbottles, rawfile)

        # Drop duplicate time variable
        ds = ds.drop('timeS')
        
        # Write dataset to new netcdf file
        ds.to_netcdf(nc_full, 'a', 'NETCDF4')

        print("File %s.nc saved successfully\n" % nc_file)

    except Exception as e:
        
        print(str(e) + 'for ' + nc_file)

        # Delete file if creation started
        if rootcheck:
            rootgrp.close()
            os.remove(nc_full)
            print("""Error saving file %s \nFile deleted\n""" % nc_file)

    return None


        
### Code Input ###

fieldtrip = 'RS2019'
##fieldtrip = 'KISSME2017'
##transect = 'T9'
transect = None



##def main():
if 1==1:
    
    print('Executing main')
    
    # Specify work locations
    workDir = os.path.abspath(os.curdir)[:-12]
    
    # Define work directories and find CTD profile data
    rawDir = os.path.join(workDir, 'data', 'field', fieldtrip, 'raw_data', 'profiling')
    saveDir = rawDir.replace('raw_data', 'processed_data')
    # Check is save dir exists and make if not
    if not os.path.isdir(saveDir):
        os.makedirs(saveDir)

    # Read field trip CTD logbook
    logFileName = fieldtrip + '_CTD_Logbook.xlsx' 
    logbook = os.path.join(workDir, 'fieldwork', fieldtrip, logFileName)
    df_check = pd.read_excel(logbook)

    ## Loop through all folders in raw data dir including rawDir
    if transect is not None:
        folds = [os.path.join(rawDir, transect)]

    else:
        x = os.listdir(rawDir)
        folds = [x[0] for x in os.walk(rawDir)]

    for fol in folds:
        
        # Find all CNV files
        rsFileSearch = os.path.join(fol, '*.cnv')
        cnv_files = glob.glob(rsFileSearch)


        ## Loop through file list
        for file in cnv_files:

            # Check and rename cnv file
            f_str = os.path.split(file)
            if f_str[-1][0:4] == 'CTD_':
                f_new = os.path.join(f_str[0], f_str[-1][4:])
                os.rename(file, f_new)
                # Rename bl, hdr, hex, XMLCON files as well
                os.rename((file[0:-3] + 'bl'), (f_new[0:-3] + 'bl'))
                os.rename((file[0:-3] + 'hdr'), (f_new[0:-3] + 'hdr'))
                os.rename((file[0:-3] + 'hex'), (f_new[0:-3] + 'hex'))
                os.rename((file[0:-3] + 'XMLCON'), (f_new[0:-3] + 'XMLCON'))
                file = f_new

            # Check and rename cnv file
            if fieldtrip == 'RS2019':
                if os.path.split(file)[-1][0] != 'L':
                    f_str = os.path.split(file)
                    if f_str[-1][0] != 'T':
                        tloc = f_str[-1].find('T')
                        f_new = os.path.join(f_str[0], (f_str[-1][tloc:tloc+3] + \
                                                        f_str[-1][0:tloc] + f_str[-1][tloc+3:]))
                        os.rename(file, f_new)
                        # Rename bl, hdr, hex, XMLCON files as well
                        os.rename((file[0:-3] + 'bl'), (f_new[0:-3] + 'bl'))
                        os.rename((file[0:-3] + 'hdr'), (f_new[0:-3] + 'hdr'))
                        os.rename((file[0:-3] + 'hex'), (f_new[0:-3] + 'hex'))
                        os.rename((file[0:-3] + 'XMLCON'), (f_new[0:-3] + 'XMLCON'))
                        file = f_new

            # Set nc filename from cnv name
            saveFin = fol.replace('raw_data', 'processed_data')
            # Check is save dir exists and make if not
            if not os.path.isdir(saveFin):
                os.makedirs(saveFin)
            nc_file, nc_full = extra.set_NC_filename(file, saveFin)

            # Read CNV file
            sbe = seabird.CnvFile(file, edit=False)

            # Read CNV header
            header = read_cnv_header(file)

            # Create CTD timeseries
            time_series = create_CTD_timeseries(sbe, header)

            # Get units list
            units = get_SBE_units(sbe)

            # Check and ammend bottle position array
            sbe = check_bottle_array(sbe)

            # Check if valid filename in logbook
##            df_check_temp = df_check[df_check['Filename'].str.match(nc_file[:-3], na=False)]

            # Build xarray of variables
            ds = xr.Dataset()

            # Insert datenum variable first
            V = xr.DataArray(time_series, dims=('depSM',), name='datenum',\
                    attrs = {'units':'days since 1970-01-01', 'longname': 'Date number'},\
                    coords = {'depSM':sbe.array[:,sbe.names.index('depSM')].tolist()})

            ds.update({'datenum':V})

            # Convert SBE to xarray dataset
            ds = sbe_to_dataset(ds, sbe, units)


            ## Create NetCDF
            create_ctd_nc(nc_full, saveFin, ds, sbe, file, fieldtrip)

##        return None
    

##if __name__ == '__main__':
##    main()

