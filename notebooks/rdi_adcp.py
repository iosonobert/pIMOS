# coding: utf-8

"""
Process RDI adcp data using dolfyn (lkilcher.github.io/dolfyn/)
"""

#from dolfyn.adp import api
from dolfyn.io import rdi
#from dolfyn.adp import rotate

# AZ commenting these out because my cython is not working. 
# from transform import Transform, rdi_xyz_enu

import xarray as xr
import os
from datetime import datetime, timedelta
from collections import OrderedDict
import numpy as np
from scipy.interpolate import interp1d

from matplotlib.dates import num2date, date2num

import rdi_adcp_utils
import pdb

deg2rad = np.pi / 180.
rad2deg = 180./np.pi

class rdi_adcp(object):
    """
    Wrapper object for parsing an RDI file with the "dolfyn" package from github

    Attributes:
    --------
        ds   - xr dataset with the pertinent data
        _data - adcp data object from dolfyn

    Methods:
    --------
        to_netcdf : save xr.Dataset to netcdf

    """

    
    def __init__(self, adcpfile, group=None, rotate=False, \
        mapbins=False, config_update = {}):
        """
        Inputs:
        -----
                adcpfile - RDI binary or a netcdf file
                rotate - [False] option to rotate
                config_update - dictionary with keys to replace in the ADCP config file

        """

        self.adcpfile = adcpfile

        self.rotate = rotate
        self.mapbins = mapbins


        try:
            # Load data directly from a netcdf file
            self.ds = xr.open_dataset(self.adcpfile, group=group, decode_coords=False)

        except:
            # Load the raw binary data
            self._data = rdi.read_rdi(self.adcpfile) # This would be Levi Kilcher's dat object

            # Override settings in the configuration dictionary
            for kk in list(config_update.keys()):
                self._data.config[kk] = config_update[kk]

            self.positive = self._data.config['orientation']
                
            # Rotate raw data
            if rotate:
                print('Rotating...')
                #self._rotate_velocity()
                #self._rotate_velocity_uhdas(self._data._u, \
                self.u, self.v, self.w, self.errvel,\
                self.u_inst, self.v_inst, self.w_inst = \
                    self._rotate_velocity(self._data.vel.astype('float64'), \
                        self._data.range,\
                        self._data.orient.raw.heading.astype('float64'),\
                        self._data.orient.raw.pitch.astype('float64'),\
                        self._data.orient.raw.roll.astype('float64'),\
                        float(self._data.config['beam_angle']),\
                        self._data.config['beam_pattern'],\
                        self._data.config['orientation'],\
                        self._data.config['coord_sys'])

            self.time = num2date(self._data.mpltime)

            self.ds, self.encoding = self._to_xray(self._data)

    def to_netcdf(self, ncfile, mode='w', group=None):
        """
        Save to a netcdf file
        """
        print('Saving data to file %s...'%ncfile)

        if 'encoding' in list(self.__dict__.keys()):
            encoding=self.encoding
        else:
            encoding = None

        self.ds.to_netcdf(ncfile,mode=mode, group=group, encoding=encoding)
        print('Done.')

    #######
    # Private methods
    #######
    def _calc_depth(self, orientation=None, P=None):
        """
        Calculate the bin depths as a 2D array using the instrument pressure
        """
        ds = self.ds

        tilt = self._calc_tilt(ds['pitch'].values,\
               ds['roll'].values)

        if orientation is None:
            orientation = ds.orientation

        # Calculate the depth using the TILT and pressure
        z = ds.distance.values
        dz = np.abs(z[1] - z[0])
        zhat = z[:,np.newaxis] * np.cos(tilt*deg2rad)

        # Depth of the instrument
        if P is None:
            P = ds['pressure'].values


        if ds.orientation == 'down':
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


    def _calc_tilt(self, pitch, roll):
        """
        Calculate instrument tilt
        """

        for i in np.arange(0, 10):
            print('MAKE THIS A PROPERTY')

        # See the IMOS wiki for recomendations on this:
        # https://github.com/aodn/imos-toolbox/wiki/QCProcedures#adcp-tilt-test---imostiltvelocitysetqc---compulsory
        #TILT = acos(sqrt(1 - sin(ROLL)^2 - sin(PITCH)^2))

        TILT = np.arccos(np.sqrt(1 - np.sin(roll*deg2rad)**2\
                - np.sin(pitch*deg2rad)**2)) * rad2deg

        return TILT

    def _rotate_velocity(self, beamvel, distance,\
        heading_deg, pitch_deg, roll_deg,\
        beam_angle, beam_pattern, orientation, coord_sys):
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
        if self.mapbins:
            beamvelnew = binmap(beamvel, distance, pitch_deg, roll_deg)
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

    def _to_xray(self, data):
        """
        Create a dictionary of the output variables 
        and convert to an xr Dataset
        """
        #data = self._data
        print('Converting to xarray')

        # dimensions
        adcp_dims = {
            'distance':data.range.squeeze(),
            'time':self.time,
            'beam':list(range(1,5)),
        }

        if not data.config.coord_sys == 'beam':
            raise(Exception)

        # Variables
        #      - Shift the beam to the end. I find it nicer this way so that when you pull out a single variable you have a 2D matrix - easier to plot. 

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
        }


        if self.rotate:
            adcp_vars_rotation = {
                'u':
                    {'data':self.u,
                    'attrs':{'long_name':'Eastward water velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'v':
                    {'data':self.v,
                    'attrs':{'long_name':'Northward water velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'w':
                    {'data':self.w,
                    'attrs':{'long_name':'Vertical water velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'uinst':
                    {'data':self.u_inst,
                    'attrs':{'long_name':'Instrument x velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'vinst':
                    {'data':self.v_inst,
                    'attrs':{'long_name':'Instrument y velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'winst':
                    {'data':self.w_inst,
                    'attrs':{'long_name':'Instrument upward velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                'errvel':
                    {'data':self.errvel,
                    'attrs':{'long_name':'Error velocity',
                'coordinates':'time distance',
                        'units':'m/s'},
                    'dims':('distance','time')
                    },
                    }

            print(' Adding rotation variables')
            adcp_vars.update(adcp_vars_rotation)

        var_names = list(adcp_vars.keys())

        # Create an output dataset
        # Global attrs
        attrs = data.config

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

        return ds, encoding

