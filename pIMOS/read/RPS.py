"""
Parse the RPS observations into a suitable NetCDF format
"""

import os
import operator
import xarray as xr
import numpy as np 
import pandas as pd
from datetime import datetime, timedelta


RPSdefault = {
   'SeaWaterTemp':'temperature',\
   'SeaWaterTemp1':'temperature',\
   'Salinity':'salinity',\
   'Pressure':'pressure',\
   'WaterDepth':'waterlevel',\
   'WaterDepth1':'waterlevel',\
   'SensorDepth':'waterlevel',\
   'CurEastComp':'east_vel',\
   'CurNorthComp':'north_vel',\
   'water_u':'east_vel',\
   'water_v': 'north_vel',\
   'water_w': 'up_vel',\
   'CurSpd': 'current_speed',\
   'CurDirn':'current_direction'}

RPScoords = {
   'Height':'nominal_instrument_height_asb',
   'Height1':'nominal_instrument_height_asb',
   'DepthHeight': 'nominal_site_depth',
   'Longitude':'nominal_longitude',
   'Latitude':'nominal_latitude',\
}


current_vars = {
   'SeaWaterTemp':'temperature',\
   'SeaWaterTemp1':'temperature',\
   'CurEastComp':'water_u',\
   'CurNorthComp':'water_v'}

alt_vars = {
   'SeaWaterTemp':'temperature',\
   'SeaWaterTemp1':'temperature',\
   'CurSpd':'water_vel',\
   'CurDirn':'water_dir'}


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
    if parse_times:
        qg = da.attrs['quality_group']
        Q1 = ds['Qual1'].values.astype(np.double)
        qgen = 'quality_group'
        # Matt's advice was to use flag = 1 or 2
        mask = calculate_qc_mask_old(Q1, qg, 1, 2) 

    # Calculate the flag (new format)
    else:
        mask, qgen = calculate_qc_mask_new(ds, da, verbose=verbose)

    if np.any(mask) & verbose & (qgen is not None):
        print('\t\tBad values found in ' + da.name + '... masking with ' + qgen + '!')

    # return the mask
    return mask, qgen

def add_rps_variable(ds, dsr, time, varname, newname, verbose=True, parse_times=True):
    qgen = None
    if ('Qual' not in varname) & ('quality' not in varname):
        if dsr[varname].size > 1:
            V = convert_rps_variable(dsr, time, varname, newname, verbose)
            mask, qgen = flag_rps_variable(dsr, V, parse_times, verbose)
            ds.update({newname: V})
            if qgen is not None:
                qc_var = newname + '_qc'
                qc_val = np.zeros(len(ds.time.values))
                qc_val[mask] = 1
                ds[qc_var] = xr.DataArray(qc_val, dims=('time',), name=qc_var)
                ds[newname].attrs['qc_variable'] = qc_var
    return ds, qgen


# Basic read function
def process_rps_file(filename, RPSvars=RPSdefault, parse_times=True, verbose=True, convert_all=True):
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
    dsr = xr.open_dataset(filename, decode_times=~parse_times)

    # Check for time variable
    if 'Time' not in dsr.coords.keys():
        if 'Time' not in dsr.data_vars.keys():
            raise(ValueError('No ''Time'' coordinate in file ' + filename))

    # Convert time (for old format)
    if parse_times:
        time = convert_rps_time(dsr['Time'])
    else:
        time = dsr['Time'].values

    # Create an output dataset
    ds = xr.Dataset(attrs=dsr.attrs, coords={'time':time})

    # Convert the data variables and insert them into the output object
    non_vars = []
    for varname in dsr.data_vars.keys():
        if varname != 'Time':
            print(varname)
            if varname in RPSvars.keys():
                newname = RPSvars[varname]
                ds, qgen = add_rps_variable(ds, dsr, time, varname, newname, verbose=verbose, parse_times=parse_times)
            else:
                newname = varname
                if convert_all:
                    print(dsr[varname])
                    ds, qgen = add_rps_variable(ds, dsr, time, varname, newname, verbose=verbose, parse_times=parse_times)
                    if qgen is None:
                        non_vars.append(varname)
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
    subdict['timestep'] = (nc['Time'][-1] - nc['Time'][0]).values
    subdict['end'] = nc['Time'][-1].values

    # Add location
    subdict['location'] = nc.attrs['location']
    return subdict




#######################
# Metnet functions (crap)
#######################

def load_rtcp_data(fullpath, keep='Current'):
    df = pd.read_csv(fullpath, skiprows=2, parse_dates=True, index_col=0)
    if 'Time' in df.columns:
        df = df.drop(columns='Time')
    # Drop any columns that dont contain 'current'
    df = df.loc[:, df.columns.str.contains(keep)]
    df.index = pd.to_datetime(df.index)
    return df.sort_index()

def get_RTCP_distance(df):
    dcol = df.columns.to_list()
    dcolxx = dcol[0].find('Bin') + 3
    bins = [int(d[dcolxx:dcolxx+2]) for d in dcol]
    return 4.25 + 2*np.array(bins)

def rtcp_da(df, distance, name=None): 
    return xr.DataArray(df.values, coords=[df.index, distance], dims=['time', 'distance'], name=name)

def rtcp_read_halfvar(file, varname):
    da = load_rtcp_data(file)
    dist = get_RTCP_distance(da)
    return rtcp_da(da, dist, varname)

def compile_raw_rtcp(files, dir=[0,1], spd=[2,3], temp=3, pos=4, surface=5):
    da_dir1 = rtcp_read_halfvar(files[dir[0]], 'current_dir')
    da_dir2 = rtcp_read_halfvar(files[dir[1]], 'current_dir')

    da_spd1 = rtcp_read_halfvar(files[spd[0]], 'current_spd')
    da_spd2 = rtcp_read_halfvar(files[spd[1]], 'current_spd')

    da_dir = xr.concat([da_dir1, da_dir2], dim='distance')
    da_dir.attrs['units'] = 'degree N'

    da_spd = xr.concat([da_spd1, da_spd2], dim='distance')
    da_spd.attrs['units'] = 'm s$^{-1}$'

    ds = xr.Dataset({'direction': da_dir, 'speed': da_spd})
    
    # Add seperate temp data
    df_temp = load_rtcp_data(files[temp], keep='Temp')
    ds['temperature'] = xr.DataArray(df_temp.values.flatten(), coords=[df_temp.index], dims=['time'], name='temperature')

    # Add the position data
    df_pos = load_rtcp_data(files[pos], keep='deg')
    latitude = xr.DataArray(df_pos['Buoy_Latitude/deg'], coords=[df_pos.index], dims=['time'], name='latitude')
    longitude = xr.DataArray(df_pos['Buoy_Longitude/deg'], coords=[df_pos.index], dims=['time'], name='longitude')

    latitude[latitude == 0.0] = np.nan
    longitude[longitude == 0.0] = np.nan

    ds['latitude'] = latitude.interp_like(ds['speed'])
    ds['longitude'] = longitude.interp_like(ds['speed'])

    # Add the surface CM-04
    df_surf = load_rtcp_data(files[surface], keep='Surface')
    # only add variables where index overlaps with ds
    df_surf = df_surf.loc[df_surf.index.intersection(df_temp.index)]
    ds['current_dir_surface'] = xr.DataArray(df_surf['Total_CurrentDirection/Surface_deg_Mean'],\
                                     coords=[df_surf.index], dims=['time'], name='current_dir_surface')    
    ds['current_speed_surface'] = xr.DataArray(df_surf['Total_CurrentSpeed/Surface_mps_Mean'],\
                                     coords=[df_surf.index], dims=['time'], name='current_speed_surface')     
    return ds


def read_raw_wms(file, cols=[0,1,2,3,8,9,10], colnames=None, skiprows=2,\
                 latbounds=[-13.85,-13.80], lonbounds=[123.27,123.29]):

    if colnames is None:
        colnames = ['time', 'latitude', 'longitude', 'up_vel', 'direction', 'speed', 'temperature']

    df = pd.read_csv(file, skiprows=skiprows)
    df = df.iloc[:,cols]

    df.columns = colnames

    # Remove some bad data by bounds
    ix = ((df['latitude'] > latbounds[0]) & (df['latitude'] < latbounds[1])) &\
         ((df['longitude'] > lonbounds[0]) & (df['longitude'] < lonbounds[1]))
    df = df.loc[ix]

    # Reindex to time & convert to numpy datetime64
    df.set_index('time', inplace=True)
    df.index = pd.to_datetime(df.index)

    # Convert to dataset
    ds = df.to_xarray()

    return ds