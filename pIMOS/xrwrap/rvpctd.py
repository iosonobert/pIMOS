"""
Research vessel profiling CTD wrapper.

Standalone wrapper (not inheriting from RVLCTD) that enforces a
"time-only coordinate" convention where Depth/latitude/longitude are stored
as variables rather than coordinates.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import num2date, date2num
import scipy.signal as signal
import pdb
import datetime
import os

import pIMOS.xrwrap.pimoswrap as pimoswrap
import pIMOS.read.SEABIRD_CTD as read_sbdctd

font = {'weight': 'normal', 'size': 12}
matplotlib.rc('font', **font)


def from_netcdf(infile):
    """
    Pass straight to the class-level from_netcdf loader.
    """
    rr, ds = RVPCTD.from_netcdf(infile)
    rr.add_variable_attributes()
    rr.enforce_time_only_coordinate()
    return rr, ds


def from_cnvfile(filename, depth_name='depSM'):
    """
    Spin up RVPCTD from Seabird CNV file.
    """
    ds = read_sbdctd.read(filename, depth_name=depth_name, depth_to_coord=False)

    ds.attrs['raw_file_name'] = os.path.split(filename)[1]
    ds.attrs['raw_file_directory'] = os.path.split(filename)[0]

    rr = RVPCTD(ds)

    # Toward process level 1 metadata/variable conventions.
    rr.add_variable_attributes()
    rr.enforce_time_only_coordinate()

    return rr, ds


##########################
# Actual xarray wrap #####
##########################
class RVPCTD(pimoswrap.pimoswrap):

    _default_attrs = pimoswrap.pimoswrap._default_attrs.copy()
    _default_attrs['is_profile_data'] = 1

    class_attrs = {
        'title': 'Measured data from a research vessel profiling CTD',
        'source': 'pIMOS',
    }

    def __init__(self, ds):

        print('Initialising accessor.')
        self.ds = ds  # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        self.enforce_these_attrs(self.class_attrs)

        self.enforce_time_only_coordinate()

    def add_variable_attributes(self):
        """
        Adds attributes to variables and associates QC flags to variables.
        """

        if 'depSM' in self.ds:
            self.ds = self.ds.rename(name_dict={'depSM': 'Depth'})

        if 'prDM' in self.ds:
            self.ds = self.ds.rename(name_dict={'prDM': 'Pressure'})

        if 't090C' in self.ds:
            self.ds = self.ds.rename(name_dict={'t090C': 'Temperature'})

        if 't190C' in self.ds:
            self.ds = self.ds.rename(name_dict={'t190C': 'Temperature2'})

        if 'sbeox0MmperL' in self.ds:
            self.ds = self.ds.rename(name_dict={'sbeox0MmperL': 'Oxygen'})

        if 'sbeox1MmperL' in self.ds:
            self.ds = self.ds.rename(name_dict={'sbeox1MmperL': 'Oxygen2'})

        if 'c0Sperm' in self.ds:
            self.ds = self.ds.rename(name_dict={'c0Sperm': 'Conductivity'})

        if 'c1Sperm' in self.ds:
            self.ds = self.ds.rename(name_dict={'c1Sperm': 'Conductivity2'})

        if 'density00' in self.ds:
            self.ds = self.ds.rename(name_dict={'density00': 'Density'})

        if 'flECO_AFL' in self.ds:
            self.ds = self.ds.rename(name_dict={'flECO_AFL': 'Fluorescence'})

        if 'par' in self.ds:
            self.ds = self.ds.rename(name_dict={'par': 'PAR'})

        if 'sal00' in self.ds:
            self.ds = self.ds.rename(name_dict={'sal00': 'Salinity'})

        if 'sal11' in self.ds:
            self.ds = self.ds.rename(name_dict={'sal11': 'Salinity2'})

        if 'CStarTr0' in self.ds:
            self.ds = self.ds.rename(name_dict={'CStarTr0': 'Beam_transmission'})

        if 'CStarAt0' in self.ds:
            self.ds = self.ds.rename(name_dict={'CStarAt0': 'Beam_attenuation'})

        if 'altM' in self.ds:
            self.ds = self.ds.rename(name_dict={'altM': 'Altimeter'})

        if 'accM' in self.ds:
            self.ds = self.ds.rename(name_dict={'accM': 'Acceleration'})

        if 'timeS' in self.ds and 'time' not in self.ds.coords:
            self.ds = self.ds.assign_coords({'time': self.ds.timeS})
            self.ds = self.ds.drop_vars('timeS')

        if 'TEMP' in self.ds:
            self.ds = self.ds.rename(name_dict={'TEMP': 'Temperature'})

        if 'CNDC' in self.ds:
            self.ds = self.ds.rename(name_dict={'CNDC': 'Conductivity'})

        if 'PRES' in self.ds:
            self.ds = self.ds.rename(name_dict={'PRES': 'Pressure'})

        if 'TEMP2' in self.ds:
            self.ds = self.ds.rename(name_dict={'TEMP2': 'Temperature2'})

        if 'CNDC2' in self.ds:
            self.ds = self.ds.rename(name_dict={'CNDC2': 'Conductivity2'})

        if 'sbeox0Mm_L' in self.ds:
            self.ds = self.ds.rename(name_dict={'sbeox0Mm_L': 'Oxygen'})

        if 'sbeox1Mm_L' in self.ds:
            self.ds = self.ds.rename(name_dict={'sbeox1Mm_L': 'Oxygen2'})

        if 'DEPTH' in self.ds:
            self.ds = self.ds.rename(name_dict={'DEPTH': 'Depth'})

        if 'density' in self.ds:
            self.ds = self.ds.rename(name_dict={'density': 'Density'})

        ds = self.ds

        for v in ['Temperature', 'Temperature2']:
            if v in ds.data_vars:
                ds[v].attrs['long_name'] = 'seawater_temperature'
                ds[v].attrs['standard_name'] = 'seawater_temperature'
                ds[v].attrs['units'] = 'deg'
                ds[v].attrs['cf_compliant'] = 1
                if 'qc_variable' not in ds[v].attrs:
                    self.associate_qc_flag(v, v)

        if 'Pressure' in ds.data_vars:
            ds['Pressure'].attrs['standard_name'] = 'sea_water_pressure'
            ds['Pressure'].attrs['long_name'] = '"Sea water pressure" is the pressure that exists in the medium of sea water. It includes the pressure due to overlying sea water, sea ice, air and any other medium that may be present.'
            ds['Pressure'].attrs['units'] = 'dbar'
            ds['Pressure'].attrs['cf_compliant'] = 1
            if 'qc_variable' not in ds['Pressure'].attrs:
                self.associate_qc_flag('Pressure', 'Pressure')

        for v in ['Conductivity', 'Conductivity2']:
            if v in ds.data_vars:
                ds[v].attrs['long_name'] = 'sea_water_electrical_conductivity'
                ds[v].attrs['standard_name'] = 'sea_water_electrical_conductivity'
                ds[v].attrs['units'] = 'S m-1'
                ds[v].attrs['cf_compliant'] = 1
            if 'qc_variable' not in ds[v].attrs:
                self.associate_qc_flag(v, v)

        if 'Density' in ds.data_vars:
            ds['Density'].attrs['long_name'] = 'sea_water_density'
            ds['Density'].attrs['standard_name'] = 'sea_water_density'
            ds['Density'].attrs['units'] = 'kg m-3'
            ds['Density'].attrs['cf_compliant'] = 1
            if 'qc_variable' not in ds['Density'].attrs:
                self.associate_qc_flag('Density', 'Density')

        for v in ['Salinity', 'Salinity2']:
            if v in ds.data_vars:
                ds[v].attrs['long_name'] = 'sea_water_salinity'
                ds[v].attrs['standard_name'] = 'sea_water_salinity'
                ds[v].attrs['units'] = 'PSU'
                ds[v].attrs['cf_compliant'] = 1
                if 'qc_variable' not in ds[v].attrs:
                    self.associate_qc_flag(v, v)

        if 'Depth' in ds.data_vars:
            ds['Depth'].attrs['long_name'] = 'pressure_sensor_depth_below_sea_surface'
            ds['Depth'].attrs['standard_name'] = 'pressure_sensor_depth_below_sea_surface'
            ds['Depth'].attrs['units'] = 'm'
            ds['Depth'].attrs['cf_compliant'] = 1
            if 'qc_variable' not in ds['Depth'].attrs:
                self.associate_qc_flag('Depth', 'Pressure')

        for var in ds:
            if 'cf_compliant' not in ds[var].attrs:
                print('Setting {} to be non cf_compliant'.format(var))
                ds[var].attrs['cf_compliant'] = 0

        self.enforce_time_only_coordinate()

    def enforce_time_only_coordinate(self):
        """
        Make time the sole dimension coordinate as proper numpy datetime64.

        Priority for time source:
          1. 'datenum' variable (float days since 1970-01-01) → converted to datetime64
          2. 'time' variable already in ds → promoted as-is
          3. 'timeS' variable (elapsed seconds) → promoted as-is (last resort)

        After promotion, swaps the underlying dimension to 'time' so it becomes
        a dimension coordinate, then demotes all remaining non-time non-dimension
        coordinates to plain data variables.
        """
        if 'time' not in self.ds.coords:
            if 'datenum' in self.ds:
                # Convert float days-since-1970 to proper datetime64
                time_vals = pd.to_datetime(self.ds['datenum'].values, unit='D', origin='unix')
                dim = self.ds['datenum'].dims[0]
                self.ds = self.ds.assign_coords({'time': (dim, time_vals)})
                self.ds = self.ds.drop_vars('datenum')
            elif 'time' in self.ds:
                self.ds = self.ds.assign_coords({'time': self.ds['time']})
            elif 'timeS' in self.ds:
                self.ds = self.ds.assign_coords({'time': self.ds['timeS']})
                self.ds = self.ds.drop_vars('timeS')

        # Swap dims so time is the dimension coordinate (not just a non-dim coord).
        if 'time' in self.ds.coords and 'time' not in self.ds.dims:
            time_dim = self.ds['time'].dims[0]  # e.g. 'obs'
            self.ds = self.ds.swap_dims({time_dim: 'time'})

        # Demote every non-time non-dimension coordinate to a data variable.
        non_time_coords = [c for c in list(self.ds.coords)
                           if c != 'time' and c not in self.ds.dims]
        if non_time_coords:
            self.ds = self.ds.reset_coords(non_time_coords)

    def inject_lat_lon(self, nav_time, latitude, longitude, source, lat_name='external_latitude', lon_name='external_longitude'):
        """
        Interpolate external navigation latitude/longitude onto dataset time.

        Parameters
        ----------
        nav_time : array-like
            Timestamps corresponding to the navigation samples.
        latitude : array-like
            Latitude samples on nav_time.
        longitude : array-like
            Longitude samples on nav_time.
        source : str
            Required provenance string describing the nav source.
        lat_name : str, optional
            Output latitude variable name in ds (default: 'external_latitude').
        lon_name : str, optional
            Output longitude variable name in ds (default: 'external_longitude').
        """

        if source is None or str(source).strip() == '':
            raise ValueError('A non-empty "source" is required when injecting lat/lon.')

        if 'time' not in self.ds.dims:
            raise ValueError('Dataset must have time as a dimension before injecting lat/lon.')

        nav_time = pd.to_datetime(np.asarray(nav_time))
        lat = np.asarray(latitude, dtype=float)
        lon = np.asarray(longitude, dtype=float)

        if not (len(nav_time) == len(lat) == len(lon)):
            raise ValueError('nav_time, latitude, and longitude must have the same length.')

        # Keep only finite nav points with valid timestamps.
        valid = (~pd.isna(nav_time)) & np.isfinite(lat) & np.isfinite(lon)
        nav_time = nav_time[valid]
        lat = lat[valid]
        lon = lon[valid]

        if len(nav_time) < 2:
            raise ValueError('Need at least 2 valid nav points to interpolate lat/lon onto time.')

        nav_ns = nav_time.view('int64')
        sort_idx = np.argsort(nav_ns)
        nav_ns = nav_ns[sort_idx]
        lat = lat[sort_idx]
        lon = lon[sort_idx]

        # Remove duplicate timestamps so np.interp has a strictly increasing axis.
        nav_ns_unique, unique_idx = np.unique(nav_ns, return_index=True)
        lat_unique = lat[unique_idx]
        lon_unique = lon[unique_idx]

        if len(nav_ns_unique) < 2:
            raise ValueError('Need at least 2 unique nav timestamps to interpolate lat/lon.')

        target_time = pd.to_datetime(self.ds['time'].values)
        target_ns = target_time.view('int64')

        lat_interp = np.interp(target_ns, nav_ns_unique, lat_unique, left=np.nan, right=np.nan)
        lon_interp = np.interp(target_ns, nav_ns_unique, lon_unique, left=np.nan, right=np.nan)

        self.ds[lat_name] = (('time',), lat_interp)
        self.ds[lon_name] = (('time',), lon_interp)

        self.ds[lat_name].attrs['long_name'] = 'latitude'
        self.ds[lat_name].attrs['standard_name'] = 'latitude'
        self.ds[lat_name].attrs['units'] = 'degrees_north'

        self.ds[lon_name].attrs['long_name'] = 'longitude'
        self.ds[lon_name].attrs['standard_name'] = 'longitude'
        self.ds[lon_name].attrs['units'] = 'degrees_east'

        nav_comment = (
            'Injected {lat_name}/{lon_name} by linear interpolation onto ds.time '
            'from source "{source}" using {n} valid nav samples.'
        ).format(lat_name=lat_name, lon_name=lon_name, source=str(source), n=len(nav_ns_unique))

        self.add_comment(str(source), nav_comment, ds_name='ds', data_var=None)
        self.add_comment(str(source), nav_comment, ds_name='ds', data_var=lat_name)
        self.add_comment(str(source), nav_comment, ds_name='ds', data_var=lon_name)


# %%
