# coding: utf-8

"""
Process RDI adcp data using dolfyn (lkilcher.github.io/dolfyn/)
"""

#from dolfyn.adp import api
from dolfyn.io import rdi
#from dolfyn.adp import rotate

# AZ commenting these out because cython is a pain. 
# from transform import Transform, rdi_xyz_enu

import xarray as xr
import os
from datetime import datetime, timedelta
from collections import OrderedDict
import numpy as np
from scipy.interpolate import interp1d

from matplotlib.dates import num2date, date2num
import zutils.xrwrap as xrwrap 
import zutils.time as ztime

from ..utils import rdi_adcp_utils as rdi_adcp_utils

import importlib

deg2rad = np.pi / 180.
rad2deg = 180./np.pi

def from_pdo(filename, rotate=False, mapbins=False):
    """
    Read PD0 file with Dolfyn. 
    """
    
    # Load the raw binary data
    data = rdi.read_rdi(filename) # This would be Levi Kilcher's dat object
    
    print(data)
    print(data.config)
    
    # self._data = data # SHOULD NOT USE THIS AS AN ATTRIBUTE AS WON'T LOAD LATER.

    if not data.config.coord_sys == 'beam':
        raise(Exception('This is coded for beam only. Update to be more flexible.'))

    for i in np.arange(0, 10):
        print('CODED FOR BEAM ONLY.')

    # Override settings in the configuration dictionary
    # for kk in list(self.config_update.keys()):
    #     self._data.config[kk] = self.config_update[kk]

    # self.positive = self._data.config['orientation']
        
    # self.time = num2date(self._data.mpltime)
    # self.time = ztime.num2date_lk(self._data.mpltime)
    time = ztime.num2date_lk(data.mpltime)

    # dimensions
    adcp_dims = {
        'distance':data.range.squeeze(),
        'time':time,
        'beam':list(range(1,5)),
    }

    adcp_vars = {
        'beamvel':
            {'data':data.vel.transpose([1, 2, 0]),
                'attrs':{'long_name':'Beam velocity',
                    'units':'m/s'},
                'dims':('distance','time','beam')
            },
        'percent_good':
            {'data':data.signal.prcnt_gd.transpose([1, 2, 0]),
                'attrs':{'long_name':'Percentage good',
                    'units':''},
                'dims':('distance','time','beam')
            },
        'echo':
            {'data':data.signal.echo.transpose([1, 2, 0]),
                'attrs':{'long_name':'Echo intensity',
                    'units':''},
                'dims':('distance','time','beam')
            },
        'corr':
            {'data':data.signal.corr.transpose([1, 2, 0]),
                'attrs':{'long_name':'Correlation',
                    'units':''},
                'dims':('distance','time','beam')
            },
        'pressure':
            {'data':data.depth_m,
                'attrs':{'long_name':'Pressure',
                    'units':'decibars'},
                'dims':('time',)
            },
        'temperature':
            {'data':data.env.temperature_C,
                'attrs':{'long_name':'Water temperature',
                    'units':'degrees C'},
                'dims':('time',)
            },
        'heading':
            {'data':data.orient.raw.heading,
                'attrs':{'long_name':'Instrument heading',
                    'units':'degrees'},
                'dims':('time',)
            },
        'pitch':
            {'data':data.orient.raw.pitch,
                'attrs':{'long_name':'Instrument pitch',
                    'units':'degrees'},
                'dims':('time',)
            },
        'roll':
            {'data':data.orient.raw.roll,
                'attrs':{'long_name':'Instrument roll',
                    'units':'degrees'},
                'dims':('time',)
            },
        'mpltime':
            {'data':data.orient.raw.roll,
                'attrs':{'long_name':'Instrument roll',
                    'units':'degrees'},
                'dims':('time',)
            },
    }

    attrs = data.config
    var_names = list(adcp_vars.keys())
    ds = xr.Dataset(attrs=attrs)

    encoding = {}
    for var in var_names:
    #for var in ['beamvel']:
        print('Converting variable: %s...'%var)
        
        coords = OrderedDict()
        for dd in adcp_vars[var]['dims']:
            coords.update({dd:adcp_dims[dd]})
            
        V = xr.DataArray( adcp_vars[var]['data'],\
            dims=adcp_vars[var]['dims'],\
            name=var,\
            attrs = adcp_vars[var]['attrs'],\
            coords = coords
        )

        ds.update({var:V})
        
        encoding.update({var:{'zlib':True,'_FillValue':-999999.}})

    # self.ds = ds
    # self.encoding = encoding

    #############################
    ## Some extra attrs that 
    ## will go to raw file attrs
    #############################
    ds.attrs['beam_angle']     = data.config['beam_angle']
    ds.attrs['beam_pattern']   = data.config['beam_pattern']
    ds.attrs['orientation']    = data.config['orientation']
    ds.attrs['coord_sys']      = data.config['coord_sys']

    ds.attrs['raw_file_name']      = os.path.split(filename)[1]
    ds.attrs['raw_file_directory'] = os.path.split(filename)[0]
    
    rr = RDI_ADCP_PD02(ds)

    #############################
    ## Higher Process level stuff
    #############################
    rr.associate_qc_flag('beamvel', 'velocity3')

    #############################################
    # Rotate raw data, mapping if necessary #####
    #############################################
    if rotate:
        rr._calc_rotations(mapbins)
        ds = rr.ds
        
    if False: # I don't know where the Dolfyn gets this info from. And thus, I don't trust this info. 
        rr.update_attribute('nominal_instrument_orientation', data.config['orientation'])

    # self.ds, self.encoding = self._to_xray(self._data)
    for i in np.arange(0, 10):
        print('DELETED THE ENCODING FROM MATT RAYSONs CODE. QUESTION FOR HIM.')

    
    return rr, ds

def from_netcdf(infile):
    """
    Pass straight to the main xrwrap from_netcdf method.
    """
   
    classhandler = RDI_ADCP_PD02

    rr, ds = xrwrap._from_netcdf(infile, classhandler)
    
    return rr, ds

##########################
# Actual xarray wrap #####
##########################
class RDI_ADCP_PD02(xrwrap.xrwrap):

    # def __init__(self, infile, rotate=False, \
    #              driver='dolfyn',\
    #              mapbins=False, config_update = {},
    #              nens=None, attributes={}):
    #     """
    #     Inputs:
    #     -----
    #             infile - File to load. Can be a string containing the full path, or a 2 element iterable with [folder, file]. 
    #                         File type will depend on the driver. RDI binary or a netcdf file.
    #             rotate - [False] option to rotate
    #             config_update - dictionary with keys to replace in the ADCP config file

    #     """

    #     self.parse_infile(infile)

    #     self.driver_name = driver

    #     spam_loader = importlib.find_loader(driver)
    #     if spam_loader is None:
    #         raise(Exception('No module found for driver {}.'.format(driver)))

    #     self.driver = importlib.import_module(driver)

    #     self.rotate = rotate
    #     self.mapbins = mapbins
    #     self.config_update = config_update # I'm not sure why this needs to be here. A question for Matt Rayson. 

    #     if driver.lower() in ['dolfyn', 'dallsporpoise']: 
    #         self.read_pd0_dd(self.folder, self.file_, attributes=attributes, nens=nens)
    #     elif driver.lower() == 'xarray':
    #         self.load(self.folder, self.file_)
    #     else:
    #         raise(Exception('{} is not a valid driver'.format(driver)))

    #     self.update_global_attrs()
    
    # def update_global_attrs(self):
    #     """
    #     Each wrapper should overload this function
    #     """

    #     # CF Compliance
    #     print('Updating attributes function of the class.')
    #     self.update_attribute('title', 'Measured data from a TDRI ADCP read from .PD0 files')
    #     self.update_attribute('institution', 'UWA')
    #     self.update_attribute('source', 'TDRI ADCP [Workhorse, Quartermaster, or Longranger]')
    #     self.update_attribute('history', '')
    #     self.update_attribute('references', '')
    #     self.update_attribute('comment', '')

    def __init__(self, ds):
        
        print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        class_attrs = {
            'title': 'Measured data from a TDRI ADCP',
            'source': 'TDRI ADCP [Workhorse, Quartermaster, or Longranger]' # Could be more specific.
        }
        self.enforce_these_attrs(class_attrs)

    # def export(self, final=False, final_folder=None, csv=True):
    #     """
    #     Overloading the base class export function.
    #     """

    #     # to NetCDF is not working after the new dolfyn update
    #     # outname = self.outname()
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

    
    #######
    # Property methods
    #######
    @property
    def tilt(self):
        """
        Instrument tilt. 
        """

        # tilt = self._calc_tilt(self.ds['pitch'].values,\
        #        self.ds['roll'].values)
        tilt = self._calc_tilt(method='combined')

        return tilt 

    #######
    # Private methods
    #######
    def _calc_rotations(self, mapbins):

        print('Rotating...')
        
        for i in np.arange(0, 10):
            print('SHOULD DO THIS FROM self.ds NOT self._data.')

        # u, v, w, errvel,\
        # u_inst, v_inst, w_inst = \
        #     self._rotate_velocity(self._data.vel.astype('float64'), \
        #         self._data.range,\
        #         self._data.orient.raw.heading.astype('float64'),\
        #         self._data.orient.raw.pitch.astype('float64'),\
        #         self._data.orient.raw.roll.astype('float64'),\
        #         float(self._data.config['beam_angle']),\
        #         self._data.config['beam_pattern'],\
        #         self._data.config['orientation'],\
        #         self._data.config['coord_sys'],
        #         mapbins)

        beamvel = self.ds['beamvel'].values.astype('float64').transpose([2, 0, 1]) # Get back to the order it wants
        u, v, w, errvel,\
        u_inst, v_inst, w_inst = \
            self._rotate_velocity(beamvel,
                self.ds['distance'].values.astype('float64'),\
                self.ds['heading'].values.astype('float64'),\
                self.ds['pitch'].values.astype('float64'),\
                self.ds['roll'].values.astype('float64'),\
                float(self.rawattrs['beam_angle']),\
                self.rawattrs['beam_pattern'],
                self.rawattrs['orientation'], 
                self.rawattrs['coord_sys'],\
                mapbins)
                
        # dimensions
        adcp_dims = {
            'distance':self.ds['distance'].values.astype('float64'),
            'time': self.ds['time'],
            'beam':list(range(1,5)),
        }

        adcp_vars_rotation = {
                'u':
                    {'data':u,
                    'attrs':{'long_name':'Eastward water velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'v':
                    {'data':v,
                    'attrs':{'long_name':'Northward water velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'w':
                    {'data':w,
                    'attrs':{'long_name':'Vertical water velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'uinst':
                    {'data':u_inst,
                    'attrs':{'long_name':'Instrument x velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'vinst':
                    {'data':v_inst,
                    'attrs':{'long_name':'Instrument y velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'winst':
                    {'data':w_inst,
                    'attrs':{'long_name':'Instrument upward velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'errvel':
                    {'data':errvel,
                    'attrs':{'long_name':'Error velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                    }

        var_names = list(adcp_vars_rotation.keys())
        ds = xr.Dataset()

        # encoding = self.encoding
        for var in var_names:
        #for var in ['beamvel']:
            print('Converting variable: %s...'%var)
            
            coords = OrderedDict()
            for dd in adcp_vars_rotation[var]['dims']:
                coords.update({dd:adcp_dims[dd]})
                
            V = xr.DataArray( adcp_vars_rotation[var]['data'],\
                dims=adcp_vars_rotation[var]['dims'],\
                name=var,\
                attrs = adcp_vars_rotation[var]['attrs'],\
                coords = coords
            )

            ds.update({var:V})
            
            # encoding.update({var:{'zlib':True,'_FillValue':-999999.}})

        ###########################################
        ### Careful here not to overwrite the attrs
        ###########################################
        attrs = self.ds.attrs
        self.ds = xr.merge((self.ds, ds))
        self.ds.attrs = attrs
        # self.encoding = encoding

        self.associate_qc_flag('u', 'velocity')
        self.associate_qc_flag('v', 'velocity')
        self.associate_qc_flag('w', 'velocity')
        self.associate_qc_flag('uinst', 'velocity')
        self.associate_qc_flag('vinst', 'velocity')
        self.associate_qc_flag('winst', 'velocity')

    def _calc_depth(self, orientation=None, P=None):
        """
        Calculate the bin depths as a 2D array using the instrument pressure
        """

        for i in np.arange(0, 10):
            print('NEED AN INTERPOLATE PRESSURE OPTION, IN THE EVENT THAT P IS TAKEN FROM ANOTHER INST.')

        for i in np.arange(0, 10):
            print('SHOULD READ ORIENTATION FROM THE ATTRIBUTES.')

        ds = self.ds

        if orientation is None:
            raise(Exception)

        # Calculate the depth using the TILT and pressure
        z = ds.distance.values
        dz = np.abs(z[1] - z[0])
        zhat = z[:,np.newaxis] * np.cos(self.tilt*deg2rad)

        # Depth of the instrument
        if type(P) == str:
            P = self.ds[P].values
        elif type(P) in [float, int]:
            P = np.ones(self.ds.time.values.shape)*float(P)

        if orientation == 'down':
            zhat = P + zhat
            #flag_depth = zhat > max_depth - 0.15*max_depth # too conservative
            #flag_depth = zhat > max_depth - 1.5*dz
        else:
            zhat = P - zhat
            #flag_depth = zhat < 0. + 0.15*max_depth
            #flag_depth = zhat < 0. + 1.5*dz

        # Update the internal object
        zda = xr.DataArray(zhat,
            dims = ('distance','time',),\
            coords = {'time':ds.time.values,\
                    'distance':ds.distance.values},
            attrs = {\
                    'long_name':'Depth below free-surface',\
                    'units':'m',\
                    'positive':'down',\
                    },
        )

        ds.update({'zhat':zda})

    def _calc_sidelobe_trim(self, P='pressure', trim=10, trim_units='%', variables=['u', 'v', 'w', 'uinst', 'vinst', 'winst']):

        rdi_adcp_utils.sidelobe_trim(self, trim=trim, trim_units=trim_units, P=P, variables=variables)
        
    def _calc_tilt(self, method='combined'):
        """
        Calculate instrument tilt
        """
        # See the IMOS wiki for recomendations on this:
        # https://github.com/aodn/imos-toolbox/wiki/QCProcedures#adcp-tilt-test---imostiltvelocitysetqc---compulsory
        # TILT = acos(sqrt(1 - sin(ROLL)^2 - sin(PITCH)^2))

        pitch = self.ds['pitch'].values
        roll = self.ds['roll'].values

        if method.lower()=='max':
            return np.maximum(abs(roll), abs(pitch))
        elif method.lower()=='combined':
            return np.arccos(np.sqrt(1 - np.sin(roll*deg2rad)**2\
                - np.sin(pitch*deg2rad)**2)) * rad2deg
        else:
            raise(Exception("Unrecognised method."))
            

    def _rotate_velocity(self, beamvel, distance,\
        heading_deg, pitch_deg, roll_deg,\
        beam_angle, beam_pattern, orientation, coord_sys, mapbins):
        """
        Rotate from beam to compass coordinates

        Using the rotation routines from the dolyfn libary

        """

        beamvel = beamvel.transpose([1, 0, 2]) #

        if not coord_sys == 'beam':
            print('Data collected in %s coordinates - not rotating.'%(coord_sys))
            # error velocity is stored first
            return beamvel[:,0,:], beamvel[:,1,:], beamvel[:,2,:]
            #return beamvel[:,1,:], -beamvel[:,0,:], beamvel[:,2,:]

        # local copy of data for convenience
        #data = self._data

        # Map onto the correct depth-bins
        if mapbins:
            beamvelnew = rdi_adcp_utils.binmap(beamvel, distance, pitch_deg, roll_deg)
        else:
            beamvelnew = beamvel

        # Calculate the rotation matrix
        #isconvex = data.config['beam_pattern'] == 'convex'
        isconvex = beam_pattern == 'convex'
        rotmat = rdi_adcp_utils.calc_beam_rotmatrix(theta = beam_angle,convex=isconvex)
        print(rotmat)

        #data.config.update({'rotmat':rotmat})

        ## Update the beam
        #data.add_data('beam1vel',data._u[:,0,:])
        #data.add_data('beam2vel',data._u[:,1,:])
        #data.add_data('beam3vel',data._u[:,2,:])
        #data.add_data('beam4vel',data._u[:,3,:])

        # Rotate the beam to instrument coordinates
        #  this populates the variables: u_inst, v_inst
        u_inst, v_inst, w_inst, errvel = rdi_adcp_utils.beam2inst(beamvelnew, rotmat)

        # Rotate the instrument to earth coordinates
        #self.u, self.v, self.w = rotate.inst2earth(data)

        u, v, w = rdi_adcp_utils.inst2earth(u_inst, v_inst, w_inst,\
            heading_deg, pitch_deg, roll_deg,\
            orientation=orientation,\
            heading_offset=None, declination=None,\
            fixed_orientation=False)

        #return u_inst, v_inst, w_inst
        return u, v, w, errvel, u_inst, v_inst, w_inst

    # def _to_xray(self, data):
    #     """
    #     Create a dictionary of the output variables 
    #     and convert to an xr Dataset
    #     """
    #     #data = self._data
    #     print('Converting to xarray')

    #     # Variables
    #     #      - Shift the beam to the end. I find it nicer this way so that when you pull out a single variable you have a 2D matrix - easier to plot. 

    #     if self.rotate:
            
    #         print(' Adding rotation variables')
    #         adcp_vars.update(adcp_vars_rotation)

    #     # Create an output dataset
    #     # Global attrs
        

    #     return ds, encoding

