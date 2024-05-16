"""
Created on Tue Apr 02 2024

@author: WE
"""

# import os 
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

########### Would need to add this to requirements ############
import pyproj as pp

import pIMOS.xrwrap.pimoswrap as xrwrap
import pIMOS.utils.UWA_archive_utils as ai 

try:
    from d2spike.despike_GN import qc0_Flags
    from d2spike.despike import D2spikearray
except:
    print('d2spike not available')
    pass

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

    # If time is not monotonic, sort it
    if not np.all(np.diff(ds.time) > np.timedelta64(0)):
        ds = ds.sortby('time')

    rr = DRIFTER(ds)

    return rr, ds


def from_netcdf(infile):
    """
    Pass straight to the main xrwrap from_netcdf method.
    """
    classhandler = DRIFTER
    rr, ds = xrwrap._from_netcdf(infile, classhandler)
    return rr, ds


def fit_spline(x, y, xnew, s=0.0):
    spl = UnivariateSpline(x, y, s=s)
    return spl(xnew)


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
        e, n = calc_east_north(self.ds[lon_var], self.ds[lat_var], proj=proj, zone=zone, ellps=ellps, south=south, **ppkwargs)
        self.ds['easting'] = xr.DataArray(e, dims='time')
        self.ds['northing'] = xr.DataArray(n, dims='time')
        return self


    # Function to calculate distance between each point in an array of lat/lon
    def calc_distance(self, lon_var='longitude', lat_var='latitude', ellps='WGS84', **ppkwargs):
        dist = calc_geodistance(self.ds[lon_var], self.ds[lat_var], ellps=ellps, **ppkwargs)
        self.ds['distance'] = xr.DataArray(dist, dims='time')
        return self
    

    def calc_raw_speed(self, dist_var='distance', time_var='time', qc_var=None):
        good_flag = qcvar(qc_var=qc_var, ds=self.ds)

        speed_c = calc_speed(self.ds[dist_var][good_flag], self.ds[time_var][good_flag])
        speed_a = np.full(len(self.ds[dist_var]), np.nan)
        speed_a[good_flag] = speed_c

        self.ds['raw_speed'] = xr.DataArray(speed_a, dims=time_var, attrs={'units':'m s$^{-1}$',\
                                                                         'comments':'First order fowards in time gradient ($\\delta x / \\delta t$) for QC'})
        return self
    

    def UnphysicalQC(self, var_name, flag_name, thresh=3.0, flag_before=False, flag_after=False):
        'This is mostly done in the despike algorithm, but left here for other situations.'
        qc_unphysical(self, var_name=var_name, flag_name=flag_name, thresh=thresh, flag_before=flag_before, flag_after=flag_after)
        return self
    

    def calc_vel(self, qc_var=None):
        east_vel, north_vel = calc_cart_vel(self.ds.easting, self.ds.northing, self.ds.time, qc_var=qc_var, ds=self.ds)
        self.ds['east_vel'] = xr.DataArray(east_vel, dims='time', name='eastward_sea_water_velocity',\
                                           attrs={'units':'m s$^{-1}$', 'comments':'Second order centred in time gradient ($\\delta x / \\delta t$) (np.gradient)'})
        
        self.ds['north_vel'] = xr.DataArray(north_vel, dims='time', name='northward_sea_water_velocity',\
                                           attrs={'units':'m s$^{-1}$', 'comments':'Second order centred in time gradient ($\\delta x / \\delta t$) (np.gradient)'})
        return self
    

    def despike_drifter(self, window=3, spike_dist=400, qc0_val=3.0,\
                        reinstate=0.2, timecut=None, qc_var=None):
        pos_flag = despike_drifter_position(self.ds, window=window, spike_dist=spike_dist, qc0_val=qc0_val,\
                                            reinstate=reinstate, timecut=timecut, qc_var=qc_var)

        # Calc good speed and add the flags in 
        self.associate_qc_flag('latitude', 'position')
        self.update_qc_flag_logical('qc_position', 'time', ~pos_flag, 1)
        self.update_qc_flag_logical('qc_position', 'time', self.ds['qc_position'] != 1, 0)
        self.associate_qc_flag('longitude', 'position')
        if 'easting' in self.ds:
            self.associate_qc_flag('easting', 'position')
            self.associate_qc_flag('northing', 'position')
        self.calc_vel(qc_var='qc_position')
        return self
    

    def plot_drifter(self, size=20, qc_var='qc_position', cmap_var=None, utm=False, cmap='viridis', attrs=None):
        return plot_drifter_map(self, size=size, qc_var=qc_var, cmap_var=cmap_var, utm=utm, cmap=cmap, attrs=attrs)
    


##############################################################################################


def qcvar(qc_var=None, ds=None, length=None):
    if qc_var is None:
        if length is not None:
            good_flag = np.full(length, True)
        elif ds is not None:
            good_flag = np.full(len(ds.time), True)
        else:
            raise ValueError('Either provide a length or a dataset.')
    elif isinstance(qc_var, str):
        good_flag = ds[qc_var]==0
    else:
        good_flag = qc_var
    return good_flag


def calc_east_north(lon, lat, proj='utm', zone=51, ellps='WGS84', south=True, **ppkwargs):
    p = pp.Proj(proj=proj, zone=zone, ellps=ellps, south=south, **ppkwargs)
    return p(lon, lat)


def calc_geodistance(longitude, latitude, ellps='WGS84', **ppkwargs):
    p = pp.Geod(ellps=ellps, **ppkwargs)
    dist = np.full(len(longitude), np.nan)
    dist[1:] = p.inv(longitude[:-1], latitude[:-1], longitude[1:], latitude[1:])[2]
    return dist


def calc_speed(distance, time):   
    speed = np.full(len(distance), np.nan)
    speed[1:] = distance[1:] / (np.diff(time) / np.timedelta64(1, 's'))
    return speed
    

def calc_cart_vel(easting, northing, time, qc_var=None, ds=None):
    good_flag = qcvar(qc_var, ds)

    # Calculate the east-west and north-south velocity
    east_vel = np.full(len(time), np.nan)
    east_vel[good_flag] = np.gradient(easting[good_flag],\
                            ((time - time[0])[good_flag] / np.timedelta64(1, 's')))

    north_vel = np.full(len(time), np.nan)
    north_vel[good_flag] = np.gradient(northing[good_flag],\
                            ((time - time[0])[good_flag] / np.timedelta64(1, 's')))
    return east_vel, north_vel


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


def despike_drifter_position(ds, window=3, spike_dist=400, qc0_val=3.0,\
                    reinstate=0.2, timecut=None, qc_var=None):
    '''
    Despike a drifter that is an xarray dataset with variables time, easting, northing.
    '''
        
    time = (ds.time.values - ds.time.values[0]) / np.timedelta64(1, 's')
    east = ds.easting.values.copy()
    north = ds.northing.values.copy()

    # Flag time issues
    if timecut is not None:
        fx = np.where(np.diff(ds.time.values) < timecut)[0]
        east[fx] = np.nan
        north[fx] = np.nan

    # Calculate velocity with spikes included (except bad GPS data)
    east_vel, north_vel = calc_cart_vel(east, north, ds.time.values, qc_var=qc_var, ds=ds)
    east_vel = xr.DataArray(east_vel, coords={'time':ds.time.values})
    north_vel = xr.DataArray(north_vel, coords={'time':ds.time.values})
    north_vel.floatda.qc0_flags(val=qc0_val)
    east_vel.floatda.qc0_flags(val=qc0_val)

    north_background = north_vel.rolling(time=window, center=True, min_periods=window).mean()
    north_hf = north_vel - north_background

    east_background = east_vel.rolling(time=window, center=True, min_periods=window).mean()
    east_hf = east_vel - east_background

    north_hf, _ = north_hf.floatda.despike_gn23()
    east_hf, _ = east_hf.floatda.despike_gn23()

    # Fit splines to new data
    nanx = np.isnan(north_hf) | np.isnan(east_hf)
    if np.sum(~nanx) > 0:
        north_fit = fit_spline(time[~nanx], north[~nanx], time)
        east_fit = fit_spline(time[~nanx], east[~nanx], time)

        # Spike dist from line
        spike_dist_from_line = np.sqrt((north_fit - north)**2 + (east_fit - east)**2)
        east_vel, north_vel = calc_cart_vel(east_fit, north_fit, ds.time, qc_var=None, ds=ds)
    
        pos_flag = (spike_dist_from_line < spike_dist) &\
                    (np.abs(north_vel - north_background.interpolate_na('time')) < reinstate) &\
                    (np.abs(east_vel - east_background.interpolate_na('time')) < reinstate)
    else:
        pos_flag = np.ones_like(time, dtype=bool)

    return pos_flag



#### Plot utils ######################################################################################

def get_drifter_aspect(rr):
    lon_range = rr.ds.longitude.max() - rr.ds.longitude.min()
    lat_range = rr.ds.latitude.max() - rr.ds.latitude.min()
    return lon_range / lat_range


def plot_drifter_map(rr, size=20, qc_var='qc_position', cmap_var=None, utm=False, cmap='viridis', attrs=None):

    if utm:
        pos_x = 'easting'
        pos_y = 'northing'
    else:
        pos_x = 'longitude'
        pos_y = 'latitude'

    if cmap_var is None:
        # Try find speed
        if 'speed' in rr.ds:
            cmap_var = rr.ds['speed']
        else:
            cmap_var = np.sqrt(rr.ds.east_vel**2 + rr.ds.north_vel**2)
    
    fig_aspect = get_drifter_aspect(rr)
    fig, ax = plt.subplots(1, 1, figsize=(size, size/fig_aspect))

    sc = ax.scatter(rr.get_qaqc_var(pos_x),
                    rr.get_qaqc_var(pos_y),
                    c=cmap_var,
                    s=1,
                    cmap=cmap)
    plt.colorbar(sc, ax=ax, label='Speed [m s$^{-1}$]', pad=0.01)
    ax.set_aspect('equal')

    ax.scatter(rr.ds[pos_x][rr.ds['qc_position']==1],\
               rr.ds[pos_y][rr.ds['qc_position']==1],\
               s=50, c='r', marker='x')
    plt.grid()

    if attrs is None:
        attrs = ai.nonempty_attrs(rr)
    title= ' | '.join([str(attrs[i]) for i in attrs])
    title += ' | Flags=' + str(int(np.sum(rr.ds[qc_var]))) + '/' + str(len(rr.ds[qc_var]))
    plt.title(title)

    return fig, ax