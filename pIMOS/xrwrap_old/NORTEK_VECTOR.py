# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:17:42 2018

@author: 20270917
"""
#%%
import matplotlib.pyplot as plt

import pandas
import numpy as np
import xarray as xr
import matplotlib
import datetime
import os 
import pdb

import importlib 

from zutils.xrwrap import xrwrap
import zutils.stats as zstats
import zutils.file as zfile

font = {'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

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

class NORTEK_VECTOR(xrwrap):
    
    folder = ''
    file = ''
    so = ''
    eo = ''

    first_good = None
    last_good = None
    
    def __init__(self, folder=None, file_=None, attributes={}, nens=None, dat=None, driver='dolfyn'):
        """
        folder - folder where file is located
        file_  - file to read
        
        ** OPTIONAL 
            - nens - [only used if the vec driver is selected]. Number of ensembles to read 
                    - None for whole file
                    - int for an integer number of files
                    - or [int, int] for [start_ens, end_ens] 
            - dat - add the lkilcher.dolfyn dat file, if the dat driver is chosen
        """

        spam_loader = importlib.find_loader(driver)
        if spam_loader is None:
            raise(Exception('No module found for driver {}.'.format(driver)))

        self.driver = importlib.import_module(driver)
                
        if driver.lower() in ['dolfyn', 'dallsporpoise']: 
            self.read_vec_dd(folder, file_, attributes=attributes, nens=nens)
        elif driver.lower() == 'xarray':
            self.load(folder, file_)
            # self.__outname = file_ # Keep this name.
        else:
            raise(Exception('{} is not a valid driver'.format(driver)))

        self.update_global_attrs()

    def update_global_attrs(self):
        """
        Each wrapper should overload this function
        """

        # CF Compliance
        print('Updating attributes function of the class.')
        self.update_attribute('title', 'Measured data from a Nortek Vector ADV read from .VEC files')
        self.update_attribute('institution', 'UWA')
        self.update_attribute('source', 'Nortek Vector ADV')
        self.update_attribute('history', '')
        self.update_attribute('references', '')
        self.update_attribute('comment', '')

    def export(self, final=False, final_folder=None, csv=True):
        """
        Overloading the base class export function.
        """

        # to NetCDF is not working after the new dolfyn update
        outname = self.file_
        folder = self.folder

        self.parse_attributes()

        if final:
            self.ds.attrs['Disclaimer'] = self.disclaimer
            outname = outname + 'finalised'
            if not final_folder is None:
                folder = final_folder

        self.ds.close() # Force close

        if csv:
            self.ds.to_dataframe().to_csv('{folder}//{file_}.csv'.format(folder=folder, file_=outname))

        nc_file = '{folder}//{file_}.nc'.format(folder=folder, file_=outname)
        self.ds.to_netcdf(path=nc_file)

        return self.ds

    def read_vec_dd(self, folder, file_, nens=None, attributes={}):
        """
        Read vector with either Dolfyn of Dall's porpoise. 
        """

        self.folder = folder
        self.file_ = file_
        
        if not type(self.file_) == list:
            self.file_ = [self.file_]
                
        for i in np.arange(0, len(self.file_)):

            print(self.file_)
            print(i)

            ds_, dat_ = self.read_vec_dd_single(i, nens=nens)

            # New version of dolfyn outputs timezone aware pandas timestamp objects 
            # this breaks everything so we must convert to np.datetime64
            # Do this in this statment for backward copatibility with old versions of Dolfyn
            if type(ds_.time.values[0]) is pandas._libs.tslibs.timestamps.Timestamp:
                print('Must convert date to timezone naive format')
                time = ds_.time.values
                new_time = np.array([i.to_datetime64() for i in time])
                ds_ = ds_.assign_coords({'time': new_time})
                print('time zone converted')

            # if not ds_.time.values[0].tzinfo is None:
            #     print('Must convert date to timezone naive format')
            #     time = ds_.time.values
            #     new_time = np.array([i.replace(tzinfo=None) for i in time])
            #     ds_ = ds_.assign_coords({'time': new_time})
            #     print('time zone converted')
                
            if i == 0:
                ds = ds_
                dat = dat_
            else:
                try:
                    dat.append(dat_) # Klicher's code is a bit iffy here
                except:
                    pass
                ds = xr.concat((ds, ds_), 'time')

        self.update_global_attrs()
        self.update_attributes_with_dict(attributes)

        coords = dat_.config.user.CoordSystem
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

        self.dat = dat
        self.ds = ds

        # Associate some QC flags
        self.associate_qc_flag('vel_dolfyn', 'velocity')
        self.associate_qc_flag('temperature', 'temperature')
        self.associate_qc_flag('pressure', 'pressure')
        self.associate_qc_flag('heading', 'compass')
        self.associate_qc_flag('pitch', 'tilt')
        self.associate_qc_flag('roll', 'tilt')

    def read_vec_dd_single(self, i, nens=None):

        # Read with either dolfyn or dallsporpoise 
        
        if not nens is None:
            print('Reading {} ensembles'.format(nens))
        else:
            print('Reading whole file')

        # should check drivers are valid based on what is actually installed. 
        print('Reading vec with {}'.format(self.driver))

        # Get away from this multiple file reading
        # dat = self.driver.read(self.fullpath(i), nens=nens)
        dat = self.driver.read(self.fullpath, nens=nens)

        print('Converting mpltime to date. Expect delays!')
        time = dat.mpltime
        time_datetime = num2date_lk(time[0]) + np.cumsum(np.array([datetime.timedelta(days=0)] + [datetime.timedelta(days=x) for x in np.diff(time)]))
        print('Converted. How were those delays?')
        
        ds_export = xr.Dataset({'vel_dolfyn': (['beam', 'time'], dat.vel), 
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
    
        ds_diag = xr.Dataset({},
                        coords={'beam': np.arange(1, 4),'time': time})

        ds = xr.merge([ds_export, ds_diag])

        return ds, dat

    def to_pto(self, date_lims=[None, None], pitch=None, roll=None, heading=None, ori='up'):

        """
        This is really inefficient. Need to get datatimes into the reading by default. Do this with timedeltas, should be fast. 
        
        SHOULD GET ORI FROM THE FILE!! 
        
        """    

        adv_object = importlib.import_module('turbo_tools.classes.adv_object') 
        
        ds_ = self.ds
        time = self.ds.datetime.values

        if not None in date_lims:
            ind1 = np.where(time>=np.datetime64(date_lims[0]))[0]
            ind2 = np.where(time<np.datetime64(date_lims[1]))[0]
            index = np.arange(ind1[0], ind2[-1]+1)
            ds_ = ds_.isel({'time':index})
        
        pto = adv_object.PointTurbulenceClass(ds_.time, 
                                ds_.vel_dolfyn[0, :], 
                                ds_.vel_dolfyn[1, :], 
                                ds_.vel_dolfyn[2, :], 
                                fs=ds_.attrs['config:fs'],
                                set_time=ds_.datetime.values)
        
        # Works for fixed instruments only
        if heading is None:
            heading = float(np.median(ds_.heading.values))
        if pitch is None:
            pitch = float(np.median(ds_.pitch.values))
        if roll is None:
            roll = float(np.median(ds_['roll'].values))
        
        pto.set_orientation(heading, pitch, roll, ori)
        
        return pto

    def clean(self):
        """
        For this, initialise a turbo_tools point turbulence class, do the necessary calcs, and return the outpiuts to self.ds. 
        """

        raise(Exception('Import from turbo_tools/notebooks required'))
        
    def turbulence_calcs(self):
        """
        For this, initialise a turbo_tools point turbulence class, do the necessary calcs, and return the outputs to self.ds. 
        """

        raise(Exception('Import from turbo_tools/notebooks required'))

    @property
    def ds_diag(self):
        """
        Simple decorator to pull only certain variables from the full dataset. 
        """

        return self.ds[['heading', 'pitch', 'roll', 'echo', 'corr', 'pressure']]


#%%
