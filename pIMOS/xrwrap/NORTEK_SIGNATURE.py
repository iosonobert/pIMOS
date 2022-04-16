# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:17:42 2018

@author: 20270917
"""
#%%
import xarray as xr
from turbo_tools.utils import time as turbo_utils_time
import matplotlib.pyplot as plt

import importlib
import datetime
import turbo_tools
import numpy as np
from turbo_tools.classes.adcp_object import TurbulenceProfilerClass
import dolfyn as dlfn

from os import listdir
import os

from zutils.time import num2date_lk as num2date_lk
import zutils.xrwrap as xrwrap

from pIMOS.utils.nortek_signature_utils import beam2inst, inst2earth

deg2rad = np.pi / 180.
rad2deg = 180./np.pi

def from_ad2cp(filename, nens=None, dat=None, hasb=0):
    """
    Initialise from raw ad2cp file. 
    """

    fullpath = filename

    if dat is None:
        if not nens is None: 
            # Before reading the full fuls, read the first and last ensembles in nens to and print the limits. 

            start, end = nens

            dat_ = dlfn.read(fullpath, nens=[start, start+2], rebuild_index=False)
            st = dat_.mpltime[0]

            dat_ = dlfn.read(fullpath, nens=[end, end+2], rebuild_index=False)
            en = dat_.mpltime[0]

            print('Ensemble range spans {} to {}'.format(num2date_lk(st), num2date_lk(en)))

            print('     {} hours in record.'.format((en-st)*24))

        dat = dlfn.read(fullpath, nens=nens, rebuild_index=False)
    else:
        pass
    
    st = dat.mpltime[0]
    en = dat.mpltime[-1]
    print('File read, {} hours in record.'.format((en-st)*24))
    
    # Some checks here. Ideally these would all be handled doen the line, but this is not the case at present. 
    if not "mpltime_b5" in dat:
        raise(Exception("4 beam files not handled."))

    b5_lag = dat.mpltime-dat.mpltime_b5
    if any(b5_lag>1/15):
        raise(Exception("Something wrong with b5 lag."))

    if not len(dat.range)==len(dat.range_b5) or not all(dat.range==dat.range_b5):
        raise(Exception("Handling range!= range_b5 is not handled."))

    if not dat.config.coord_sys=="BEAM" or not dat.config.coord_sys_b5=="BEAM":
        raise(Exception("BEAM data only handled."))


    nc = dat.vel.shape[1]
    nt = dat.vel.shape[2]

    blank1 = np.nan*np.zeros((nc, nt))
    blank3 = np.nan*np.zeros((3, nc, nt))

    blank1.shape

    print('Converting mpltime to date. Expect delays!')
    time = dat.mpltime
    time_datetime = num2date_lk(time[0]) + np.cumsum(np.array([datetime.timedelta(days=0)] + [datetime.timedelta(days=x) for x in np.diff(time)]))
    print('Converted.')
    print('     Starts: {}'.format(time_datetime[0]))
    print('     Ends: {}'.format(time_datetime[-1]))
    print('     How were those delays?')

    nc = dat.vel.shape[1]
    nt = dat.vel.shape[2]

    blank1 = np.nan*np.zeros((nc, nt))
    blank3 = np.nan*np.zeros((3, nc, nt))

    vel_dolfyn = np.concatenate((dat.vel, dat.vel_b5))
    echo = np.concatenate((dat.signal.amp, dat.signal.amp_b5))
    corr = np.concatenate((dat.signal.corr, dat.signal.corr_b5))

    ds = xr.Dataset({'vel_dolfyn': (['beam', 'height', 'time'], vel_dolfyn),
                'vel_enu': (['cartesian_axes', 'height', 'time'], blank3),
                'vel_xyz': (['cartesian_axes', 'height', 'time'], blank3),
                'four_beam_error_velocity': (['height', 'time'], blank1),
                'pressure': ('time', dat.env.press),
                'temperature': ('time', dat.env.temp),
                'speed_of_sound': ('time', dat.env.c_sound),
                'heading': (['time'], dat.orient.raw.heading),
                'pitch': (['time'], dat.orient.raw.pitch),
                'roll': (['time'], dat.orient.raw.roll),
                'echo': (['beam', 'height', 'time'], echo), 
                'corr': (['beam', 'height', 'time'], corr),
                'mpltime': ('time', dat.mpltime),
                'distance': ('height', np.arange(0, nc)),
                'cell': ('height', np.arange(0, nc))},
                coords={'beam': np.arange(1, 6), 'cartesian_axes': np.arange(1, 4), 'height': np.arange(0, nc), 'time': time_datetime})

    ds.attrs['raw_file_name']      = os.path.split(filename)[1]
    ds.attrs['raw_file_directory'] = os.path.split(filename)[0]
    
    coords = dat.config.coord_sys                
    ds.vel_dolfyn.attrs['standard_name'] = '{}_seawater_velocity'.format(coords) 
    ds.vel_dolfyn.attrs['long_name'] = 'Velocity output by Dolfyn package when reading the .vec file with no key word arguments. Coord sys was {}'.format(coords) 
    ds.vel_dolfyn.attrs['units'] = 'm/s' 

    ds.vel_xyz.attrs['standard_name'] = 'XYZZ_seawater_velocity'
    ds.vel_xyz.attrs['long_name'] = 'Instrument coordinate [XYZZ] velocity. Conventions are nortek conventions.'
    ds.vel_xyz.attrs['units'] = 'm/s' 

    ds.vel_enu.attrs['standard_name'] = 'ENU_seawater_velocity'
    ds.vel_enu.attrs['long_name'] = 'Earth coordinate [ENU] velocity. '
    ds.vel_enu.attrs['units'] = 'm/s' 

    ds.pressure.attrs['standard_name'] = 'pressure' 
    ds.pressure.attrs['long_name'] = 'pressure' 
    ds.pressure.attrs['units'] = 'dbar' 

    ds.temperature.attrs['standard_name'] = 'seawater_temperature' 
    ds.temperature.attrs['long_name'] = 'seawater_temperature' 
    ds.temperature.attrs['units'] = 'deg' 

    ds.speed_of_sound.attrs['standard_name'] = 'speed_of_sound' 
    ds.speed_of_sound.attrs['long_name'] = 'Speed of sound used by the instrument to estimate range etc.' 
    ds.speed_of_sound.attrs['units'] = 'm/s' 

    ds.distance.attrs['standard_name'] = 'distance_between_cell_and_instrument' 
    ds.distance.attrs['long_name'] = 'distance between cell and instrument, positive away from instrument' 
    ds.distance.attrs['units'] = 'm' 

    ds.height.attrs['standard_name'] = 'height_above_seabed' 
    ds.height.attrs['long_name'] = 'Height above seabed, positive upwards' 
    ds.height.attrs['units'] = 'm' 

    # Specify variables which will not be CF compliant. 
    ds['heading'].attrs['cf_compliant'] = 0
    ds['pitch'].attrs['cf_compliant'] = 0
    ds['roll'].attrs['cf_compliant'] = 0
    ds['echo'].attrs['cf_compliant'] = 0
    ds['corr'].attrs['cf_compliant'] = 0
    ds.time.attrs['cf_compliant'] = 0
    ds.beam.attrs['cf_compliant'] = 0

    attrs = {}
    attrs['config:SerialNum'] = dat.config.SerialNum
    attrs['config:SerialNum_b5'] = dat.config.SerialNum_b5
    attrs['config:fs'] = dat.config.fs
    attrs['config:blanking'] = dat.config.blanking
    attrs['config:blanking_b5'] = dat.config.blanking_b5
    attrs['config:cell_size'] = dat.config.cell_size
    attrs['config:cell_size_b5'] = dat.config.cell_size_b5
    attrs['config:coord_sys'] = dat.config.coord_sys
    attrs['config:coord_sys_b5'] = dat.config.coord_sys_b5
    attrs['config:data_desc'] = dat.config.data_desc
    attrs['config:data_desc_b5'] = dat.config.data_desc_b5
    attrs['config:model'] = dat.config.model
    attrs['config:nbeams'] = dat.config.nbeams
    attrs['config:nbeams_b5'] = dat.config.nbeams_b5
    attrs['config:ncells'] = dat.config.ncells
    attrs['config:ncells_b5'] = dat.config.ncells_b5
    attrs['config:nom_corr'] = dat.config.nom_corr
    attrs['config:nom_corr_b5'] = dat.config.nom_corr_b5
    attrs['config:power_level'] = dat.config.power_level
    attrs['config:power_level_b5'] = dat.config.power_level_b5
    attrs['config:vel_scale'] = dat.config.vel_scale
    attrs['config:vel_scale_b5'] = dat.config.vel_scale_b5

    attrs['config:TransMatrix'] = np.array2string(dat.config.TransMatrix)
    attrs['config:head:TransMatrix_howToTeadInNumpy'] = "T=T.replace('[','');T=T.replace(']','');T=np.fromstring(T, dtype=float, sep=' ').reshape((4, 4))"

    ds['distance'] = ds.cell*attrs['config:cell_size']+attrs['config:blanking']
    ds['height'] = ds.cell*attrs['config:cell_size']+attrs['config:blanking'] + hasb # This assumes upward looking. 
    
    ds['echo'] = ds.echo.astype(float)
    ds.attrs = attrs

    rr = NORTEK_SIGNATURE(ds)

    #############################
    ## Higher Process level stuff
    #############################

    # Associate some QC flags
    rr.associate_qc_flag('vel_dolfyn', 'velocity')
    rr.associate_qc_flag('vel_enu', 'velocity3')
    rr.associate_qc_flag('vel_xyz', 'velocity3')
    rr.associate_qc_flag('temperature', 'temperature')
    rr.associate_qc_flag('pressure', 'pressure')
    rr.associate_qc_flag('heading', 'compass')
    rr.associate_qc_flag('pitch', 'tilt')
    rr.associate_qc_flag('roll', 'tilt')

    return rr, rr.ds   



class NORTEK_SIGNATURE(xrwrap.xrwrap):
    
    # def __init__(self, infile, driver='dolfyn', nens=None, dat=None, hasb=0, attributes={}):
    #     """
    #     hasb is the height above seabed
    #     """
    #     self.parse_infile(infile)
        
    #     self.driver_name = driver
    #     print(driver)
        
    #     spam_loader = importlib.find_loader(driver)
    #     if spam_loader is None:
    #         raise(Exception('No module found for driver {}.'.format(driver)))

    #     self.driver = importlib.import_module(driver)

    #     if driver.lower() in ['dolfyn']: 
    #         self.read_d(dat=dat, nens=nens, hasb=0, attributes=attributes)
    #     elif driver.lower() == 'xarray':
    #         self.load(self.folder, self.file_)
    #     else:
    #         raise(Exception('{} is not a valid driver'.format(driver)))

    #     self.update_global_attrs()

    def __init__(self, ds):
        
        print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        class_attrs = {
            'title': 'Measured data from a Nortek Signature',
            'source': 'Nortek Signature ADCP' # Could be more specific.
        }
        self.enforce_these_attrs(class_attrs)
        
    
    
    # def export(self, final=False, final_folder=None, csv=True):
    #     """
    #     Overloading the base class export function.
    #     """

    #     # to NetCDF is not working after the new dolfyn update
    #     outname = self.file_
    #     folder = self.folder

    #     self.parse_attributes()

    #     if final:
    #         self.ds.attrs['Disclaimer'] = self.disclaimer
    #         outname = outname + 'finalised'
    #         if not final_folder is None:
    #             folder = final_folder

    #     self.ds.close() # Force close

    #     if csv:
    #         self.ds.to_dataframe().to_csv('{folder}//{file_}.csv'.format(folder=folder, file_=outname))

    #     nc_file = '{folder}//{file_}.nc'.format(folder=folder, file_=outname)
    #     self.ds.to_netcdf(path=nc_file)

    #     return self.ds
        
    def _calc_rotations(self):
        """
        Calculate the earth and instrument velocities. This should follow beamwise QC.   
        """

        beam_vel = self.get_qaqc_var('vel_dolfyn')[0:4, :, :].values
        XYZZ = beam2inst(beam_vel, self.T)
        self.ds.vel_xyz.values = XYZZ[0:3, :, :]

        hh = (self.ds['heading'])
        pp = (self.ds['pitch'])
        rr = (self.ds['roll'])

        ENU = inst2earth(self.ds.vel_xyz.values, hh, pp, rr,\
                orientation='up',
                heading_offset=None, declination=None,
                fixed_orientation=True)

        self.ds.vel_enu.values = ENU
        
        print('Done.')
        
    @property
    def T(self):
        """
        Alias for TransMatrix
        """
        
        return self.TransMatrix
        
    @property
    def TransMatrix(self):
        """
        Get the 4 BEAM BEAM to INST Transformation Matrix
        """
        
        T  = self.get_raw_file_attributes()['config:TransMatrix']
        T = T.replace('[','');T=T.replace(']','');T=np.fromstring(T, dtype=float, sep=' ').reshape((4, 4))
        
        return T

    def _calc_tilt(self, method='combined'):
        """
        Calculate instrument tilt
        """

        pitch = self.ds['pitch'].values
        roll = self.ds['roll'].values

        if method.lower()=='max':
            return np.maximum(abs(roll), abs(pitch))
        elif method.lower()=='combined':
            return np.arccos(np.sqrt(1 - np.sin(roll*deg2rad)**2\
                - np.sin(pitch*deg2rad)**2)) * rad2deg
        else:
            raise(Exception("Unrecognised method."))
  

# signature = NORTEK_SIGNATURE(fullpath, dat=signature.dat, nens=nens)
