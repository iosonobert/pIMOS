# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:17:42 2018

@author: 20270917
"""
#%%
# import matplotlib.pyplot as plt

import pandas
import numpy as np
import xarray as xr
import matplotlib
import datetime
import os 
# import pdb
import importlib 

import pIMOS.xrwrap.pimoswrap as pimoswrap

# import zutils.xrwrap as xrwrap
# import zutils.file as zfile
# import zutils.stats as zstats
# import pIMOS.xrwrap.xrwrap as xrwrap
# import pIMOS.utils.file as zfile

font = {'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

def from_vec_v0(filename, nens=None, driver='dolfyn', debug=False):
    """
    Read vector with either Dolfyn v0 or Dall's porpoise. 
    """
    
    if type(filename) == list:
        raise(Exception("Reading muyltiple files no longer supported"))
            
    
    # This was used prior to 2023
    ds, dat = read_vec_dd_single(filename, driver=driver, nens=nens, debug=debug)

    # This was put in in 2023
    # ds, dat = read_vec_dd_single_DV1(filename, nens=nens, debug=debug)

    # New version of dolfyn outputs timezone aware pandas timestamp objects 
    # this breaks everything so we must convert to np.datetime64
    # Do this in this statment for backward copatibility with old versions of Dolfyn
    if type(ds.time.values[0]) is pandas._libs.tslibs.timestamps.Timestamp:
        print('Must convert date to timezone naive format')
        time = ds.time.values
        new_time = np.array([i.to_datetime64() for i in time])
        ds = ds.assign_coords({'time': new_time})
        print('time zone converted')

    # if not ds_.time.values[0].tzinfo is None:
    #     print('Must convert date to timezone naive format')
    #     time = ds_.time.values
    #     new_time = np.array([i.replace(tzinfo=None) for i in time])
    #     ds_ = ds_.assign_coords({'time': new_time})
    #     print('time zone converted')
            
    coords = dat.config.user.CoordSystem
    ds.vel_dolfyn.attrs['standard_name'] = '{}_seawater_velocity'.format(coords) 
    ds.vel_dolfyn.attrs['long_name'] = 'Velocity output by Dolfyn package when reading the .vec file with no key word arguments. Coord sys was {}'.format(coords) 
    ds.vel_dolfyn.attrs['units'] = 'm/s' 

    ds.pressure.attrs['standard_name'] = 'pressure' 
    ds.pressure.attrs['long_name'] = 'pressure' 
    ds.pressure.attrs['units'] = 'dbar' 

    ds.temperature.attrs['standard_name'] = 'seawater_temperature' 
    ds.temperature.attrs['long_name'] = 'seawater_temperature' 
    ds.temperature.attrs['units'] = 'deg' 
    
    ds.speed_of_sound.attrs['standard_name'] = 'speed_of_sound' 
    ds.speed_of_sound.attrs['long_name'] = 'Speed of sound used by the instrument to estimate range etc.' 
    ds.speed_of_sound.attrs['units'] = 'm/s' 

    # Specify variables which will not be CF compliant. 
    ds['heading'].attrs['cf_compliant'] = 0
    ds['pitch'].attrs['cf_compliant'] = 0
    ds['roll'].attrs['cf_compliant'] = 0
    ds['echo'].attrs['cf_compliant'] = 0
    ds['corr'].attrs['cf_compliant'] = 0
    ds.time.attrs['cf_compliant'] = 0
    ds.beam.attrs['cf_compliant'] = 0

    attrs = {}
    attrs['config:fs'] = dat.config.fs

    # THESE DON'T SEEM TO ALWAYS BE THERE. 
    # attrs['config:checkdata:First_samp'] = dat.config.checkdata.First_samp
    # attrs['config:checkdata:Samples'] = dat.config.checkdata.Samples
    # attrs['config:data_header:Corr1'] = dat.config.data_header.Corr1
    # attrs['config:data_header:Corr2'] = dat.config.data_header.Corr2
    # attrs['config:data_header:Corr3'] = dat.config.data_header.Corr3
    # attrs['config:data_header:NRecords'] = dat.config.data_header.NRecords
    # attrs['config:data_header:Noise1'] = dat.config.data_header.Noise1
    # attrs['config:data_header:Noise2'] = dat.config.data_header.Noise2
    # attrs['config:data_header:Noise3'] = dat.config.data_header.Noise3
    # attrs['config:data_header:time'] = dat.config.data_header.time
    # THESE DON'T SEEM TO ALWAYS BE THERE. 
    
    attrs['config:hardware:FWversion'] = dat.config.hardware.FWversion
    attrs['config:hardware:HWrevision'] = dat.config.hardware.HWrevision
    attrs['config:hardware:PICversion'] = dat.config.hardware.PICversion
    attrs['config:hardware:ProLogFWver'] = dat.config.hardware.ProLogFWver
    attrs['config:hardware:ProLogID'] = dat.config.hardware.ProLogID
    attrs['config:hardware:config'] = dat.config.hardware.config
    attrs['config:hardware:recSize'] = dat.config.hardware.recSize
    attrs['config:hardware:serialNum'] = dat.config.hardware.serialNum
    attrs['config:hardware:status'] = dat.config.hardware.status
    attrs['config:head:NBeams'] = dat.config.head.NBeams
    attrs['config:head:TransMatrix'] = np.array2string(dat.config.head.TransMatrix)
    attrs['config:head:TransMatrix_howToTeadInNumpy'] = "T=T.replace('[','');T=T.replace(']','');T=np.fromstring(T, dtype=float, sep=' ').reshape((3, 3))"
    attrs['config:head:config'] = dat.config.head.config
    attrs['config:head:freq'] = dat.config.head.freq
    attrs['config:head:serialNum'] = dat.config.head.serialNum
    attrs['config:head:config'] = dat.config.head.config
    attrs['config:head:type'] = dat.config.head.type
    attrs['config:user:AdjSoundSpeed'] = dat.config.user.AdjSoundSpeed
    attrs['config:user:AnaInAddr'] = dat.config.user.AnaInAddr
    attrs['config:user:AnaOutScale'] = dat.config.user.AnaOutScale
    attrs['config:user:AvgInterval'] = dat.config.user.AvgInterval
    attrs['config:user:B0'] = dat.config.user.B0
    attrs['config:user:B1'] = dat.config.user.B1
    attrs['config:user:BinLength'] = dat.config.user.BinLength
    attrs['config:user:BurstMode'] = int(dat.config.user['Burst Mode'])
    attrs['config:user:Comments'] = dat.config.user.Comments
    attrs['config:user:CompassUpdRate'] = dat.config.user.CompassUpdRate
    attrs['config:user:CoordSystem'] = dat.config.user.CoordSystem
    attrs['config:user:CorrThresh'] = dat.config.user.CorrThresh
    attrs['config:user:DeployName'] = dat.config.user.DeployName
    attrs['config:user:DiagInterval'] = dat.config.user.DiagInterval
    attrs['config:user:DynPercPos'] = dat.config.user.DynPercPos
    attrs['config:user:MeasInterval'] = dat.config.user.MeasInterval
    attrs['config:user:NBeams'] = dat.config.user.NBeams
    attrs['config:user:NBeamsCellDiag'] = dat.config.user.NBeamsCellDiag
    attrs['config:user:NBins'] = dat.config.user.NBins
    attrs['config:user:NBurst'] = dat.config.user.NBurst
    attrs['config:user:NPingsDiag'] = dat.config.user.NPingsDiag
    attrs['config:user:NSamp'] = dat.config.user.NSamp
    attrs['config:user:NSampDiag'] = dat.config.user.NSampDiag
    attrs['config:user:Npings'] = dat.config.user.Npings
    attrs['config:user:ProfileTiming'] = dat.config.user['Profile Timing']
    # attrs['config:user:PwrCtrlReg'] = dat.config.user.PwrCtrlReg
    # attrs['config:user:QualConst'] = dat.config.user.QualConst
    attrs['config:user:SWVersion'] = dat.config.user.SWVersion
    attrs['config:user:SampleonSync'] = int(dat.config.user['Sample on Sync'])
    attrs['config:user:StartonSync'] = int(dat.config.user['Start on Sync'])
    attrs['config:user:WrapMode'] = dat.config.user.WrapMode
    attrs['config:user:Transmit:blankdistance'] = dat.config.user.Transmit['blank distance']
    attrs['config:user:Transmit:pulselength'] = dat.config.user.Transmit['pulse length']
    attrs['config:user:Transmit:receivelength'] = dat.config.user.Transmit['receive length']
    attrs['config:user:Transmit:time_between_bursts'] = dat.config.user.Transmit['time_between_bursts']
    attrs['config:user:Transmit:time_between_pings'] = dat.config.user.Transmit['time_between_pings']
    attrs['config:user:mode:analog_output_mode'] = int(dat.config.user.mode['analog_output_mode'])
    attrs['config:user:mode:cell_position'] = dat.config.user.mode['cell_position']
    attrs['config:user:mode:output_format'] = dat.config.user.mode['output_format']
    attrs['config:user:mode:rate'] = dat.config.user.mode['rate']
    attrs['config:user:mode:rate'] = int(dat.config.user.mode['serial_output'])
    attrs['config:user:mode:vel_scale'] = dat.config.user.mode['vel_scale']

    ds.attrs = attrs

    ds.attrs['raw_file_name']      = os.path.split(filename)[1]
    ds.attrs['raw_file_directory'] = os.path.split(filename)[0]
    
    rr = NORTEK_VECTOR(ds)
    rr.dat = dat

    # Associate some QC flags
    rr.associate_qc_flag('vel_dolfyn', 'velocity')
    rr.associate_qc_flag('temperature', 'temperature')
    rr.associate_qc_flag('pressure', 'pressure')
    rr.associate_qc_flag('heading', 'compass')
    rr.associate_qc_flag('pitch', 'tilt')
    rr.associate_qc_flag('roll', 'tilt')

    return rr, ds

def from_vec_v1(filename, nens=None, debug=False):
    """
    Read vector with Dolfyn v1. 
    """
    
    if type(filename) == list:
        raise(Exception("Reading muyltiple files no longer supported"))
            
    # This was put in in 2023
    ds, dsd = read_vec_dd_single_DV1(filename, nens=nens, debug=debug)

    coords = dsd.attrs['coord_sys']
    ds.vel_dolfyn.attrs['standard_name'] = '{}_seawater_velocity'.format(coords) 
    ds.vel_dolfyn.attrs['long_name'] = 'Velocity output by Dolfyn package when reading the .vec file with no key word arguments. Coord sys was {}'.format(coords) 
    ds.vel_dolfyn.attrs['units'] = 'm/s' 

    ds.pressure.attrs['standard_name'] = 'pressure' 
    ds.pressure.attrs['long_name'] = 'pressure' 
    ds.pressure.attrs['units'] = 'dbar' 

    ds.temperature.attrs['standard_name'] = 'seawater_temperature' 
    ds.temperature.attrs['long_name'] = 'seawater_temperature' 
    ds.temperature.attrs['units'] = 'deg' 
    
    ds.speed_of_sound.attrs['standard_name'] = 'speed_of_sound' 
    ds.speed_of_sound.attrs['long_name'] = 'Speed of sound used by the instrument to estimate range etc.' 
    ds.speed_of_sound.attrs['units'] = 'm/s' 

    # Specify variables which will not be CF compliant. 
    ds['heading'].attrs['cf_compliant'] = 0
    ds['pitch'].attrs['cf_compliant'] = 0
    ds['roll'].attrs['cf_compliant'] = 0
    ds['echo'].attrs['cf_compliant'] = 0
    ds['corr'].attrs['cf_compliant'] = 0
    ds.time.attrs['cf_compliant'] = 0
    ds.beam.attrs['cf_compliant'] = 0

    attrs = {}
    attrs['config:fs'] = dsd.attrs['fs']
    

    ds.attrs = attrs
    
    rr = NORTEK_VECTOR(ds)
    
    ds.attrs['raw_file_name']      = os.path.split(filename)[1]
    ds.attrs['raw_file_directory'] = os.path.split(filename)[0]

    rr.dsd = dsd

    # Associate some QC flags
    rr.associate_qc_flag('vel_dolfyn', 'velocity')
    rr.associate_qc_flag('temperature', 'temperature')
    rr.associate_qc_flag('pressure', 'pressure')
    rr.associate_qc_flag('heading', 'compass')
    rr.associate_qc_flag('pitch', 'tilt')
    rr.associate_qc_flag('roll', 'tilt')

    return rr, ds


def read_vec_dd_single_DV1(filename, nens=None, debug=True):

    # Read with Dolfyn version 1.X.X 
    
    if not nens is None:
        print('Reading {} ensembles'.format(nens))
    else:
        print('Reading whole file')

    dolfyn = importlib.import_module('dolfyn')

    if not dolfyn.__version__[0:2] == '1.':
        raise(Exception('Valid for dolfyn version 1 only.'))

    # Get away from this multiple file reading
    # dat = self.driver.read(self.fullpath(i), nens=nens)
    dsd = dolfyn.read(filename, nens=nens, debug=debug)

    ds = xr.Dataset({'vel_dolfyn': (['beam', 'time'], dsd.vel.values), 
                        'pressure': ('time', dsd.pressure.values),
                        'temperature': ('time', dsd.temp.values),
                        'speed_of_sound': ('time', dsd.c_sound.values),
                        'heading': (['time'], dsd.heading.values), 
                        'pitch': (['time'], dsd.pitch.values), 
                        'roll': (['time'], dsd['roll'].values), 
                        'echo': (['beam', 'time'], dsd.amp.values), 
                        'corr': (['beam', 'time'], dsd.corr.values)},
                        coords={'beam': np.arange(1, 4), 'time': dsd.time.values})
    
    return ds, dsd

def read_vec_dd_single(filename, driver, nens=None, debug=True):

    # Read with either dolfyn or dallsporpoise 
    
    if not nens is None:
        print('Reading {} ensembles'.format(nens))
    else:
        print('Reading whole file')

    # should check drivers are valid based on what is actually installed. 
    print('Reading vec with {}'.format(driver))

    spam_loader = importlib.find_loader(driver)
    if spam_loader is None:
        raise(Exception('No module found for driver {}.'.format(driver)))

    driver = importlib.import_module(driver)

    # Get away from this multiple file reading
    # dat = self.driver.read(self.fullpath(i), nens=nens)
    dat = driver.read(filename, nens=nens, debug=debug)

    print('Converting mpltime to date. Expect delays!')
    time = dat.mpltime
    time_datetime = num2date_lk(time[0]) + np.cumsum(np.array([datetime.timedelta(days=0)] + [datetime.timedelta(days=x) for x in np.diff(time)]))
    print('Converted. How were those delays?')
    
    ds = xr.Dataset({'vel_dolfyn': (['beam', 'time'], dat.vel), 
                    'pressure': ('time', dat.env.pressure),
                    'temperature': ('time', dat.env.temp),
                    'speed_of_sound': ('time', dat.env.c_sound),
                    'heading': (['time'], dat.orient.raw.heading), 
                    'pitch': (['time'], dat.orient.raw.pitch), 
                    'roll': (['time'], dat.orient.raw.roll), 
                    'echo': (['beam', 'time'], dat.signal.amp), 
                    'corr': (['beam', 'time'], dat.signal.corr),
                    'datetime': ('time', time_datetime)},
                    coords={'beam': np.arange(1, 4), 'time': time})

    return ds, dat
    
class NORTEK_VECTOR(pimoswrap.pimoswrap):


    class_attrs = {
            'title': 'Measured data from a Nortek Vector',
            'source': 'pIMOS' 
        }

    def __init__(self, ds, verbose=True):
        self.verbose = verbose

        if self.verbose:
            print('Initialising accessor.')  
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)
        self.update_attributes_with_dict(class_attrs)
        self.enforce_these_attrs(class_attrs)

    def to_pto(self, pitch=None, roll=None, heading=None, ori='up', date_lims=[None, None], **kwargs):
        """
        Return a turbo_lance point turbulence object for turbulence calcs etc. 

        NOTE: turbo_lance point turbulence objects currently only support fixed instruments

        Inputs:
            pitch: [float]. Constant pitch to use for the rotations. Default is None, in which case the record average will be used.   
            roll:  [float]. Constant roll to use for the rotations. Default is None, in which case the record average will be used. 
            roll:  [float]. Constant roll to use for the rotations. Default is None, in which case the record average will be used. 
            ori:   [str]{'up', 'down}. Orientation of the vector [the body, not the head] to use for the rotations. Default is 'up'.  

        Additional inputs for further processing.
             phase_unwrap: [bool] whether or not to include phase unwrapping
             despike:      [bool] whether or not to include despiking
             calc_mean:    [bool] whether or not to calculate mean quantities
        """    

        adv_object = importlib.import_module('turbo_lance.classes.adv_object') 
        
        ds_ = self.ds
        time = self.ds.datetime.values

        if not None in date_lims:
            ind1 = np.where(time>=np.datetime64(date_lims[0]))[0]
            ind2 = np.where(time<np.datetime64(date_lims[1]))[0]
            index = np.arange(ind1[0], ind2[-1]+1)
            ds_ = ds_.isel({'time':index})
        
        fs = self.get_raw_file_attributes()['config:fs']

        pto = adv_object.PointTurbulenceClass(ds_.time.values, 
                                ds_.vel_dolfyn[0, :].values, 
                                ds_.vel_dolfyn[1, :].values, 
                                ds_.vel_dolfyn[2, :].values, 
                                fs=fs,
                                set_time=ds_.datetime.values)
        
        # Works for fixed instruments only
        if heading is None:
            heading = float(np.median(ds_.heading.values))
        if pitch is None:
            pitch = float(np.median(ds_.pitch.values))
        if roll is None:
            roll = float(np.median(ds_['roll'].values))
        
        pto.set_orientation(heading, pitch, roll, ori)


        ########################
        # ADDITIONAL PROCESING #
        ########################
        phase_unwrap = kwargs.pop('phase_unwrap', False)
        despike      = kwargs.pop('despike', False)
        calc_mean    = kwargs.pop('calc_mean', False)

        if phase_unwrap:
            print('Running unwrap')
            pto.clean_unwrap(block_length_seconds=30)
        if despike:
            print('Running despike')
            pto.clean_despike(block_length_seconds=120)
        if calc_mean:
            print('Running mean calcs')
            pto.calculate_mean_enu()
        
        print('PTO complete')
        return pto

    def clean(self):
        """
        For this, initialise a turbo_lance point turbulence class, do the necessary calcs, and return the outpiuts to self.ds. 
        """

        raise(Exception('Import from turbo_lance/notebooks required'))
        
    def turbulence_calcs(self):
        """
        For this, initialise a turbo_lance point turbulence class, do the necessary calcs, and return the outputs to self.ds. 
        """

        raise(Exception('Import from turbo_lance/notebooks required'))

    @property
    def ds_diag(self):
        """
        Simple decorator to pull only certain variables from the full dataset. 
        """

        return self.ds[['heading', 'pitch', 'roll', 'echo', 'corr', 'pressure']]

#########################
# MISC FUNCTIONS ########
#########################
def num2date_lk(mpltime):
    """
    Reimplementation of datetime num2date to stop date errors. That occur continually when updating dolfyn or datetime. 
    
    This was copied from the Levi Kilcher's Dolfyn package, hence the lk suffix. 
    """

    for ii in np.arange(0, 5):
        print('SHOULD IMPORT ZUTILS.TIME.num2date_lk')

    if isinstance(mpltime, np.ndarray):
        try:
            n = len(mpltime)
        except TypeError:
            pass
        else:
            out = np.empty(n, dtype='O')
            for idx, val in enumerate(mpltime.flat):
                out[idx] = num2date_lk(val)
            out.shape = mpltime.shape
            return out
        
    if np.isnan(mpltime):
        return None
    return datetime.datetime.fromordinal(int(mpltime)) + datetime.timedelta(days=mpltime % 1)

#%%
