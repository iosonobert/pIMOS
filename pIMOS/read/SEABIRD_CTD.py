"""
    Tools for working with Seabird CTD files.

    Some code from UHDAS

    Some of my own...

    M.Rayson
    UWA
    Apr 2016
"""
#from __future__ import division
#from future.builtins import range
#from future.builtins import object
#from future.builtins import PY3
# We are using the native str() builtin, so don't import it from future.

import re
import datetime
import numpy as np

from datetime import datetime, timedelta
import xarray as xray
import os

import pdb

from pIMOS.utils import othertime
from pIMOS.utils import seabird_utils
import xarray as xr 

def read(fullpath, edit=False):
    
    # Read CNV file
    sbe = seabird_utils.CnvFile(fullpath, edit=edit)

    # Read CNV header
    header = read_cnv_header(fullpath)

    # Get units list
    units = get_SBE_units(sbe)

    # Check and ammend bottle position array
    sbe = check_bottle_array(sbe)

    # Convert SBE to xarray dataset
    ds = sbe_to_dataset(sbe, header, units)
    
    return ds

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

def sbe_to_dataset(sbe, header, units):

    '''Convert all arrays in an SBE class to arrays in a dataset, ds,
        with dimensions depSM and units specified by a list of unit strings
    '''

    # Create CTD timeseries
    time_series = create_CTD_timeseries(sbe, header)

    # Build xarray of variables
    ds = xr.Dataset()

    # Insert datenum variable first
    V = xr.DataArray(time_series, dims=('depSM',), name='datenum',\
            attrs = {'units':'days since 1970-01-01', 'longname': 'Date number'},\
            coords = {'depSM':sbe.array[:,sbe.names.index('depSM')].tolist()})

    ds.update({'datenum':V})

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
