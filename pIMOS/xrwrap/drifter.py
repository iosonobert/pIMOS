"""
Created on Tue Apr 02 2024

@author: WE
"""

# import os 
import numpy as np
import pandas as pd
import xarray as xr

########### Would need to add this to requirements ############
import pyproj as pp

import pIMOS.xrwrap.pimoswrap as xrwrap

class_attrs = {
            'title': 'GPS position data from a drifter',
            'source': 'pIMOS',
            'process_level': 1,
            'is_profile_data': 1,
            'nominal_latitude':np.nan,
            'nominal_longitude':np.nan
        }


def from_csv(file, timevar=None, keepcols=None, **kwargs):
    """
    Read a drifter csv file into an xarray dataset

    Parameters
    ----------
    file : str
        Path to the csv file
    timevar : str, optional
        Name of the time variable, by default None
    keepcols : list, optional
        List of columns to keep, by default None
    **kwargs : dict
        Additional arguments to pass to pd.read_csv
        
    Returns
    -------
    xarray.Dataset
        xarray dataset
    """
    # Read csv file into a pandas dataframe
    df = pd.read_csv(file, **kwargs)
    # df = pd.read_csv(file, headerlines=headerlines, parse_dates=True, index_col=index_col, encoding=encoding)

    # Convert to xarray
    if keepcols:
        # Check if list of strings or indices
        if isinstance(keepcols[0], str):
            ds = df.loc[:, keepcols].to_xarray()
        elif isinstance(keepcols[0], int):
            ds = df.iloc[:, keepcols].to_xarray()
    else:
        ds = df.to_xarray()

    # Rename coord to Time
    if timevar:
        ds = ds.rename({timevar: 'time'})
    else:
        # Rename the index to Time if not already
        if 'time' not in ds.coords:
            try:
                ds = ds.rename({'index': 'time'})
            except:
                print('No time variable found, and index is not a time variable')
                pass
    ds = ds.drop_duplicates('time')

    # Check time and reverse if necessary
    if ds['time'][0] > ds['time'][-1]:
        ds = ds.reindex(Time=ds['time'][::-1].values)

    rr = DRIFTER(ds)

    return rr, ds


def from_netcdf(infile):
    """
    Pass straight to the main xrwrap from_netcdf method.
    """
    classhandler = DRIFTER
    rr, ds = xrwrap._from_netcdf(infile, classhandler)
    return rr, ds


##########################
# Actual xarray wrap #####
##########################
class DRIFTER(xrwrap.pimoswrap):

    def __init__(self, ds):
        
        print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        self.update_attributes_with_dict(class_attrs)
        self.drop_pIMOS_coords()

    
    def drop_pIMOS_coords(self):
        self.ds = self.ds.drop(['lat_nom', 'lon_nom', 'z_nom'])
        return self

    # Align names to CF convention
    def align_names(self, var_dict=None):
        if var_dict is not None:
            self.ds = self.ds.rename(var_dict)
        return self


    # Function to calculate Easting and Northing (makes QC easier)
    def calc_easting_northing(self, lon_var, lat_var, proj='utm', zone=51, ellps='WGS84', south=True, **ppkwargs):
        # p = pp.Proj(proj=proj, zone=zone, ellps=ellps, south=south, **ppkwargs)
        # e, n = p(self.ds[lon_var], self.ds[lat_var])
        e, n = calc_east_north(self.ds[lon_var], self.ds[lat_var], proj=proj, zone=zone, ellps=ellps, south=south, **ppkwargs)
        self.ds['Easting'] = xr.DataArray(e, dims='time')
        self.ds['Northing'] = xr.DataArray(n, dims='time')
        return self

    # Function to calculate distance between each point in an array of lat/lon
    def calc_distance(self, lon_var='longitude', lat_var='latitude', ellps='WGS84', **ppkwargs):
        p = pp.Geod(ellps=ellps, **ppkwargs)
        dist = np.full(len(self.ds[lon_var]), np.nan)
        dist[1:] = p.inv(self.ds[lon_var][:-1], self.ds[lat_var][:-1], self.ds[lon_var][1:], self.ds[lat_var][1:])[2]
        self.ds['distance'] = xr.DataArray(dist, dims='time')
        return self
    
    def calc_raw_speed(self, dist_var='distance', time_var='time'):
        speed = np.full(len(self.ds[dist_var]), np.nan)
        speed[1:] = self.ds[dist_var][1:] / (np.diff(self.ds[time_var]) / np.timedelta64(1, 's'))
        self.ds['raw_speed'] = xr.DataArray(speed, dims=time_var, attrs={'units':'m s$^{-1}$',\
                                                                         'comments':'First order fowards in time gradient ($\\delta x / \\delta t$) for QC'})
        return self
    
    def UnphysicalQC(self, var_name, flag_name, thresh=3.0, flag_before=False, flag_after=False):
        qc_unphysical(self, var_name=var_name, flag_name=flag_name, thresh=thresh, flag_before=flag_before, flag_after=flag_after)
        return self
    

    # def flag_unphysical(self, var, threshold, **kwargs):
    #     """
    #     Flag unphysical values in a variable

    #     Parameters
    #     ----------
    #     var : str
    #         Variable to flag
    #     threshold : float
    #         Threshold value for flagging
    #     **kwargs : dict
    #         Additional arguments to pass to the flagging function

    #     Returns
    #     -------
    #     xarray.Dataset
    #         xarray dataset
    #     """
    #     self.ds[var] = self.ds[var].where(self.ds[var] < threshold, np.nan)
    #     return self



def calc_east_north(lon, lat, proj='utm', zone=51, ellps='WGS84', south=True, **ppkwargs):
    p = pp.Proj(proj=proj, zone=zone, ellps=ellps, south=south, **ppkwargs)
    return p(lon, lat)

def qc_unphysical(rr, var_name, flag_name, thresh=3.0, flag_before=False, flag_after=False):
    """Quality control unphysical values in a variable."""
    ix = rr.ds[var_name] > thresh
    if flag_before:
        # Also flag the points before
        ix = ix | np.roll(ix, -1)
    if flag_after:
        # Also flag the points after
        ix = ix | np.roll(ix, 1)

    # Update the flag
    rr.update_qc_flag_logical(flag_name, 'time', ix, 1)
    return rr

# def drifter_vel