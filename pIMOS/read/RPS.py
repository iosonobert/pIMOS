"""
Parse the RPS observations into a suitable NetCDF format
"""

import os
import operator
import xarray as xr
import numpy as np 
from datetime import datetime, timedelta


RPSdefault = {
   'SeaWaterTemp':'temperature',\
   'SeaWaterTemp1':'temperature',\
   'Salinity':'salinity',\
   'Pressure':'pressure',\
   'WaterDepth':'waterlevel',\
   'SensorDepth':'waterlevel',\
   'CurEastComp':'water_u',\
   'CurNorthComp':'water_v'}

RPScoords = {
   'Height':'nominal_instrument_height_asb',
   'Height1':'nominal_instrument_height_asb',
   'DepthHeight': 'nominal_site_depth',
   'Longitude':'nominal_longitude',
   'Latitude':'nominal_latitude',\
}


############
# Functions
############

def scm(s):
    return s.lstrip().rstrip().replace(' - ','_').replace(' ','_').replace('-','_')

def lscm(s):
    return s.lower().lstrip().rstrip().replace(' - ','_').replace(' ','_').replace('-','_')

def get_arr_dict(dict, var):
    return np.array([dict[dc][var] for dc in dict])


def process_depth_attr(nc):
    depthstr = 'Depth: %f %s %s'%\
        (nc.DepthHeight, nc.DepthHeight.attrs['units'], nc.DepthHeight.attrs['height_ref'])

    return depthstr


def get_station_name(nc, stat_attr='location'):
    ### Not really sure about this one, replaced ctr with 0 for now
    site = nc.attrs[stat_attr]
    site = lscm(site)
    long_site = '%s_%04d'% (site,0)
    return site, long_site


def convert_rps_time(time):
    """
    Convert the RPS time (DataArray) into a datetime object.
    Intended for old format RPS files.

    Returns numpy datetime64 array in UTC time zone.
    """
    tzone = time.attrs['timezone']
    dtzone = timedelta(hours=float(tzone))

    tint = time.values[:,0].astype(int)
    days= tint.astype(str)
    times = [datetime.strptime(xx, '%Y%m%d')-dtzone for xx in days]

    msec = time.values[:,1].tolist()
    dt = [timedelta(milliseconds=float(mm)) for mm in msec]

    return np.array([t+Dt for t,Dt in zip(times,dt)]).astype('<M8[ns]') # datetime64 object


def calculate_qc_mask_old(Q1, qg, lower, upper):
    shift = 2**(2*qg)
    idx = Q1 < 0.0
    hibit = np.zeros_like(Q1)
    Q1[idx] = Q1[idx] + 2**31
    hibit[idx] = 1
    qual = np.mod( np.floor(Q1/shift), 4)
    if qg == 15:
        qual[idx] = qual[idx]+2.
    mask = operator.or_(qual < lower, qual > upper)
    return mask

def calculate_qc_mask_new(ds, da, verbose=True):
    fail = False
    if 'quality_variable' in da.attrs.keys():
        qgen = da.attrs['quality_variable']
        qg = ds[da.attrs['quality_variable']]
    elif 'quality_group' in da.attrs.keys():
        qgen = da.attrs['quality_group']
        qg = ds[da.attrs['quality_group']]
    elif (da.name)+'Qual' in ds.data_vars.keys():
        qg = ds[(da.name)+'Qual']
        qgen = (da.name)+'Qual'
    elif 'Qual1' in ds.data_vars.keys():
        qg = ds['Qual1']
        qgen = 'Qual1'
    else:
        fail = True
        if verbose:
            print('No QC variable found for ' + da.name)

    # Set flag to anything except 1
    if not fail:
        return qg.values != 1., qgen
    else:
        return np.ones(len(da.values)).astype('int'), None


def convert_rps_variable(ds, time, varname, newname, verbose=True):
    """
    Convert a 1-D timeseries variable.
    """
    if verbose:
        print('\tConverting variable: %s -> %s...'%(varname, newname))
    V = xr.DataArray(ds[varname].values, \
            dims=('time',),\
            name=newname,\
            attrs = ds[varname].attrs,\
            coords = {'time':time}
    )
    return V


def flag_rps_variable(ds, da, parse_times, verbose=True):
    """
    Flag a 1-D timeseries variable.
    """
    qgen = ''
    # Calculate the QC flag (old format)
    if not parse_times:
        qg = da.attrs['quality_group']
        Q1 = ds['Qual1'].values.astype(np.double)
        qgen = 'quality_group'
        # Matt's advice was to use flag = 1 or 2
        mask = calculate_qc_mask_old(Q1, qg, 1, 2) 

    # Calculate the flag (new format)
    else:
        mask, qgen = calculate_qc_mask_new(ds, da, verbose=verbose)

    if np.any(mask) & verbose:
        print('\t\tBad values found in ' + da.name + '... masking with ' + qgen + '!')

    # Null the bad values
    da[mask] = np.nan
    return da


# Basic read function
def process_rps_file(filename, RPSvars=RPSdefault, parse_times=True, verbose=True):
    '''
    Basic read file for RPS, create a new dataset with pIMOS attributes

    Parameters
    ----------
    filename : str
        The file location (full path)
    RPSvars : dict
        Dictionary to specify which RPS variables to capture and 
        how to rename them.
    parse_times : bool, optional
        A flag used to force datetime64 parsing (default is
        True). This flag indicates whether the file is in 
        the new or old RPS format.
    verbose : bool, optional
        A flag to set printing info on or off (default is on).

    Returns
    -------
    ds
        xrray DataSet with pIMOS attributes
    dsr
        xarray DataSet with original attributes (raw load)
    '''
    if verbose:
        print(72*'#')
    print('Processing file: %s'%filename)

    # Load the dataset
    dsr = xr.open_dataset(filename, decode_times=parse_times)

    # Check for time variable
    if parse_times:
        if 'Time' not in dsr.coords.keys():
            raise(ValueError('No ''Time'' coordinate in file ' + filename))
    else:
        if 'Time' not in dsr.data_vars.keys():
            if 'Time' not in dsr.coords.keys():
                raise(ValueError('No ''Time'' coordinate in file ' + filename))

    # Convert time (for old format)
    if not parse_times:
        time = convert_rps_time(dsr['Time'])
    else:
        time = dsr['Time'].values

    # Create an output dataset
    ds = xr.Dataset(attrs=dsr.attrs, coords={'time':time})

    # Convert the data variables and insert them into the output object
    non_vars = []
    for varname in dsr.data_vars.keys():
        print(varname)
        if varname in RPSvars.keys():
            newname = RPSvars[varname]
            V = convert_rps_variable(dsr, time, varname, newname, verbose)
            V = flag_rps_variable(dsr, V, parse_times, verbose)
            ds.update({newname: V})
        else:
            non_vars.append(varname)

    # Put key attributes in file
    for attr in RPScoords.keys():
        if attr in dsr.coords.keys():
            ds.attrs[RPScoords[attr]] = dsr[attr].values
        elif attr in dsr.data_vars.keys():
            ds.attrs[RPScoords[attr]] = dsr[attr].values

    # Check for failed instrument height
    if 'nominal_instrument_height_asb' not in ds.attrs.keys():
        try:
            for varname in RPSdefault.keys():
                if varname in dsr.data_vars.keys():
                    if 'height' in dsr[varname].attrs.keys():
                        ds.attrs['nominal_instrument_height_asb'] = \
                            dsr[varname].attrs['height'].astype('float')
        except:
            raise ValueError('Error looking for instrument height: ' + filename)

    # Convert site depth if necessary
    if ds.attrs['nominal_site_depth'] > 0:
        ds.attrs['nominal_site_depth'] = ds.attrs['nominal_site_depth'] * -1

    # Get the station name to attrs
    _, name = get_station_name(ds)
    ds.attrs.update({'stationname':name})
    ds.attrs.update({'original_file':filename})

    # Fix the project code (this could be moved back to the notbook level)
    if 'project' in ds.attrs:
        if 'met' in lscm(ds.attrs['project']):
            jx = ds.attrs['project'].find('J')
            ds.attrs['project'] = ds.attrs['project'][jx:jx+5]

    if verbose:
        if non_vars:
            print('Variables not added to new dataset: ' + str(non_vars))
        print('Done')
    return dsr, ds


def get_RPS_summary(nc):
    subdict = {}

    # Add data type
    subdict['type'] = nc.attrs['data_type']

    # Add position
    subdict['latitude'] = nc['Latitude'].values
    subdict['longitude'] = nc['Longitude'].values
    subdict['height'] = nc['DepthHeight'].values

    # Add record length and start time
    subdict['start'] = nc['Time'][0].values
    subdict['time'] = (nc['Time'][-1] - nc['Time'][0]).values

    # Add location
    subdict['location'] = nc.attrs['location']
    return subdict

