#from dolfyn.adp import api
from dolfyn.io import rdi
#from dolfyn.adp import rotate

# AZ commenting these out because my cython is not working. Needed only for UHDAS rotations. I.e. Not needed.
# from transform import Transform, rdi_xyz_enu

import xarray as xr
import os
from datetime import datetime, timedelta
from collections import OrderedDict
import numpy as np
from scipy.interpolate import interp1d

from matplotlib.dates import num2date, date2num

import pdb

try:
    from d2spike.despike_GN import qc0_Flags
    from d2spike.despike import D2spikearray
except:
    print('d2spike not available')
    pass

deg2rad = np.pi / 180.
rad2deg = 180./np.pi

def resample_uvw(self, dtavg, raw=False, rotate=False, othervars=[]):
    """
    Resample the velocity data by averaging and subsampling at the average interval

    Inputs:
            dtavg - average interval [seconds]
    Options:
            raw - [False] subsample raw fields
            rotate - [False] if raw, rotate as well. 
            Note that: heading/pitch/roll get averaged BEFORE transformation
            othervars - list of other variables to subsample (not u/v/w)
    """

    def resample(dset, dtavg, axis=0):
        dtstr = '%dS'%(int(dtavg))
        phi = (dset.to_pandas()).resample(dtstr, axis=axis).mean()

        # Use interpolate to fill in nan's
        try:
            return xr.DataArray(phi.interpolate(axis=axis) )
        except:
            return xr.DataArray(phi)

    # Create some global attributes
    attrs = {
            'Name':'Sub-sampled ADCP data',
            'Original File':self.adcpfile,
            'Sub-sample interval':dtavg,
    }
    
    ds = self.ds

    dsnew = xr.Dataset(attrs=attrs)

    if raw:
        dsnew.update({'beamvel' : resample( ds['beamvel'], dtavg, axis=2) })

        #####
        # Rotating after subsampling is a bad idea...
        # Leave it in here for now
        #####
        if rotate:
            # Filter the heading/pitch/roll first
            dsnew.update({'heading' : resample( ds['heading'], dtavg, axis=0) })
            dsnew.update({'pitch' : resample( ds['pitch'], dtavg, axis=0) })
            dsnew.update({'roll' : resample( ds['roll'], dtavg, axis=0) })

            u, v, w, evel, uinst, vinst, vinst = \
                self._rotate_velocity(dsnew['beamvel'], \
                    dsnew['heading'],\
                    dsnew['pitch'],\
                    dsnew['roll'],\
                    ds.attrs['beam_angle'],\
                    ds.attrs['beam_pattern'],\
                    ds.attrs['orientation'],\
                    ds.attrs['coord_sys'],)

            dsnew.update({'east_vel':u, 'north_vel':v, 'up_vel':w,'errvel':evel})
    else:
        dsnew.update({'north_vel' : resample( ds['north_vel'], dtavg, axis=1) })
        dsnew.update({'east_vel' : resample( ds['east_vel'], dtavg, axis=1) })
        dsnew.update({'up_vel' : resample( ds['up_vel'], dtavg, axis=1) })

    for vv in othervars:
            dsnew.update({vv : resample( ds[vv], dtavg, axis=0) })

    return dsnew

#######
# QA/QC methods
#######
def qaqc_pgood(self, thresh):
    """
    Mask u/v/w based on percentage good information

    Uses minimum of all four beams
    """

    raise(Exception("""
                This was from Rayson code but QC has been moved into seperate module 
                which adds comments to the file and can accomodate a more complex QC 
                flag system etc."""))

    ds = self.ds

    pgood = ds['percent_good'].min(axis=1)

    mask = pgood.values < thresh

    ds['u'].values[mask] = np.nan
    ds['v'].values[mask] = np.nan
    ds['w'].values[mask] = np.nan

    self.ds.attrs.update({'QAQC Percentage Good Threshold':thresh})

def qaqc_corr(self, thresh):
    """
    Mask u/v/w based on correlation

    Ensure that the minimum of three beams is above threshold

    Use beams 1-3 for now as these are used for the transformation

    TODO: Use the correlation to indicate which beams to use for transformation
    """
    
    raise(Exception("""
                This was from Rayson code but QC has been moved into seperate module 
                which adds comments to the file and can accomodate a more complex QC 
                flag system etc."""))

    ds = self.ds
    corr = ds['corr'][:,0:3,:].min(axis=1)

    mask = corr.values < thresh

    ds['u'].values[mask] = np.nan
    ds['v'].values[mask] = np.nan
    ds['w'].values[mask] = np.nan

    self.ds.attrs.update({'QAQC Correlation Minimum Threshold':thresh})

def qaqc_echo(self, thresh):
    """
    Mask u/v/w based on echo intensity

    Ensure that all four beams are above threshold
    """
    
    raise(Exception("""
                This was from Rayson code but QC has been moved into seperate module 
                which adds comments to the file and can accomodate a more complex QC 
                flag system etc."""))

    ds = self.ds
    echo = ds['echo'][:,:,:].min(axis=1)

    mask = echo.values < thresh

    ds['u'].values[mask] = np.nan
    ds['v'].values[mask] = np.nan
    ds['w'].values[mask] = np.nan

    self.ds.attrs.update({'QAQC Echo Intensity Minimum Threshold':thresh})

def qaqc_errvel(self, thresh):
    """
    Mask u/v/w based on error velocity

    """

    raise(Exception("""
                This was from Rayson code but QC has been moved into seperate module 
                which adds comments to the file and can accomodate a more complex QC 
                flag system etc."""))

    ds = self.ds

    evel = ds['errvel']

    mask = evel.values > thresh

    ds['east_vel'].values[mask] = np.nan
    ds['north_vel'].values[mask] = np.nan
    ds['up_vel'].values[mask] = np.nan

    self.ds.attrs.update({'QAQC Error Velocity Maximum Threshold':thresh})

def qaqc_tilt(self, cutoff_angle):
    """
    Mask u/v/w arrays when the instrument tilt goes too far
    """
    
    raise(Exception("""
                This was from Rayson code but QC has been moved into seperate module 
                which adds comments to the file and can accomodate a more complex QC 
                flag system etc."""))

    ds = self.ds

    tilt = self._calc_tilt(ds['pitch'].values,\
            ds['roll'].values)

    flag_angle = tilt > cutoff_angle

    # Create a 2d mask array
    nz = ds.distance.shape
    mask = flag_angle[np.newaxis,:].repeat(nz,0)

    ds['east_vel'].values[mask] = np.nan
    ds['north_vel'].values[mask] = np.nan
    ds['up_vel'].values[mask] = np.nan

    self.ds.attrs.update({'QAQC Maximum Tilt Angle':cutoff_angle})

def qaqc_depth(self, max_depth, orientation=None, P=None, variables=['east_vel', 'north_vel', 'up_vel']):
    """
    Mask out regions outside of the maximum depth
    """

    raise(Exception("""
                This was from Rayson code but QC has been moved into seperate module 
                which adds comments to the file and can accomodate a more complex QC 
                flag system etc."""))

    ds = self.ds
    if orientation is None:
        orientation = ds.orientation
    else:
        self.ds.attrs['orientation'] = orientation # Update this

    if 'zhat' not in list(ds.keys()):
        self._calc_depth(orientation=orientation, P=P)

    zhat = self.ds.zhat.values
    z = ds.distance.values
    dz = np.abs(z[1] - z[0])

    if orientation == 'down':
        #flag_depth = zhat > max_depth - 0.15*max_depth # too conservative
        mask = zhat > max_depth - 1.5*dz
    else:
        #flag_depth = zhat < 0. + 0.15*max_depth
        mask = zhat < 0. + 1.5*dz

    for variable in variables:
        self.ds[variable].values[mask] = np.nan
    # ds['u'].values[mask] = np.nan
    # ds['v'].values[mask] = np.nan
    # ds['w'].values[mask] = np.nan

def sidelobe_trim(self, trim=1.5, P='depth', variables=['east_vel', 'north_vel', 'up_vel'], trim_units='native'):
    """
    This actually sets these variables to nan. It does not flag them. 
    """

    sidelobe_trim_calc(self, trim=trim, P=P, trim_units=trim_units)
    for variable in variables:
        # This should be a flag update not nan
        # self.ds[variable].values[self.ds['sidelobe_blank'].values] = np.nan
        qc_var = self.ds[variable].attrs['qc_variable']
        self.update_qc_flag_logical(qc_var, 'time', self.ds['sidelobe_blank'].values, 1)

def sidelobe_trim_calc(self, trim=1.5, P='depth', trim_units='native', egvar='east_vel'):
    """
    inputs:
        P - the variable to trim by. Can be numeric or the name of a variable [e.g. pressure] with the appropriate dimensions.
               -- for upward facing instruments the 'pressure' variable is a good choice, or a nominal pressure if the pressure sensor failed.  
               -- for downward facing instruments the instrument height above the bed would be a good choice, or a nominal height if this is unknown.
        """

    # Make the blanking variables
    if not 'sidelobe_last_good' in self.ds.data_vars.keys():
        self.ds['sidelobe_last_good'] = xr.DataArray(
            data=np.zeros_like(self.ds.time),
            dims=self.ds.time.dims,
            coords=self.ds.time.coords,
            attrs=dict())
        
    if not 'sidelobe_blank' in self.ds.data_vars.keys():
        self.ds['sidelobe_blank'] = xr.DataArray(
            data=np.zeros_like(self.ds[egvar]),
            dims=self.ds[egvar].dims,
            coords=self.ds[egvar].coords,
            attrs=dict(),)    
        
    if type(P) == str:
        P = self.ds[P].values
    elif type(P) in [float, int]:
        P = np.ones(self.ds.time.values.shape)*float(P)
    else:
        raise(Exception('P can only be a variable name in the dataset, or a float.'))
        
    s = self.ds['sidelobe_blank'].shape
    self.ds['sidelobe_blank'][:] = 0
    self.ds['sidelobe_last_good'][:] = 0

    for step in np.arange(0, s[1]):
        if trim_units.lower() == 'native':
            pass
            trim_ = trim
        elif trim_units.lower() in ['percentage', 'perc', '%']:
            trim_ = P[step]*trim/100
            
        i = self.ds['distance'].values > P[step]-trim_
        self.ds['sidelobe_blank'][i, step] = 1

        ni = np.invert(i)
        if sum(ni) > 0:
            # print(self.ds['distance'].values[ni])
            self.ds['sidelobe_last_good'][step] = np.max(self.ds['distance'].values[ni])

    self.ds['sidelobe_blank'].values = self.ds['sidelobe_blank'].values != 0

    print('Trimmed sidelobe')


def despike_adcp(rr, vars_ds=['east_vel','north_vel','up_vel'], re_val=0.3, qc0_val=2.0, gf_sig=1.0, sw_thresh=0.9, verbose=False):
    for var in vars_ds:
        # Check for attribute
        if 'qc_variable' in rr.ds[var].attrs:
            qc_var = rr.ds[var].attrs['qc_variable']
            if verbose:
                print(f'Flags in {qc_var} before spiking {var} = {np.sum(rr.ds[qc_var].values)}')
        else:
            rr.associate_qc_flag(var, 'velocity')
            qc_var = 'qc_velocity'

        # Get the data
        flag_data = rr.get_qaqc_var(var)
        flag_data = flag_data.transpose('time','distance')
        orig_data = flag_data.copy()

        # Flag unphysical
        flag_data = flag_data.floatda.qc0_flags(val=qc0_val)

        # Calculate background
        data_bg = flag_data.floatda.gaussian_filter(gf_sig)

        # Subtract the background values and despike
        data_hf = flag_data - data_bg
        for ii, wd in enumerate(data_hf.T):
            data_hf[:,ii], _ = wd.floatda.despike_gn23(sw_thresh=sw_thresh, verbose=verbose)

        # Call 2D indexing reinstatement
        re_ix = np.abs(orig_data - orig_data.where(~np.isnan(data_hf)).interpolate_na('time')) < re_val
        data_hf = data_hf.floatda.reinstate_threshold((orig_data - data_bg), re_ix)
        new_data = data_hf + data_bg

        # Calculate background
        data_bg = new_data.floatda.gaussian_filter(gf_sig)

        # Subtract the background values and despike
        data_hf = new_data - data_bg
        for ii, wd in enumerate(data_hf.T):
            data_hf[:,ii], _ = wd.floatda.despike_gn23(sw_thresh=sw_thresh, verbose=verbose)

        # Call 2D indexing reinstatement
        re_ix = np.abs(orig_data - orig_data.where(~np.isnan(data_hf)).interpolate_na('time')) < re_val/2
        data_hf = data_hf.floatda.reinstate_threshold((orig_data - data_bg), re_ix)

        # Set the flag
        rr.update_qc_flag_logical('qc_velocity', 'time', np.isnan(data_hf), 1)
        if verbose:
            print(f'Flags in {qc_var} after spiking {var} = {np.sum(rr.ds[qc_var].values)}')
    return rr
    

#########
# Instrument rotation routines
#
# Adapted from dolfyn library's .../adp/rotate.py
#########

def calc_rotation_terms(theta, convex, degrees):
    if degrees:
        theta = theta * deg2rad
    if convex:
        c = 1.
    else:
        c = -1.

    a = 1 / (2. * np.sin(theta))
    b = 1 / (4. * np.cos(theta))
    d = a / (2. ** 0.5)

    return a, b, c, d

def calc_beam_rotmatrix(theta=20.0, convex=True, degrees=True):
    """Calculate the rotation matrix from beam coordinates to
    instrument head coordinates.

    Parameters
    ----------
    theta : is the angle of the heads (usually 20 or 30 degrees)

    convex : is a flag for convex or concave head configuration.

    degrees : is a flag which specifies whether theta is in degrees
        or radians (default: degrees=True)
    """
    a,b,c,d = calc_rotation_terms(theta, convex, degrees)

    return np.array([[c * a, -c * a, 0, 0],
                     [0, 0, -c * a, c * a],
                     [b, b, b, b],
                     [d, d, -d, -d]])

    # Consistent with bm2dir.m for an upward looking instrument
    #return np.array([[-c * a, c * a, 0, 0],
    #                 [0, 0, -c * a, c * a],
    #                 [-b, -b, -b, -b],
    #                 [-d, -d, d, d]])

def beam2inst(beamvel, rotmat):
    """Rotate velocities from beam to instrument coordinates.
    """
    #if hasattr(adcpo.config, 'rotmat'):
    #    rotmat = adcpo.config.rotmat
    #else:
    #    rotmat = calc_beam_rotmatrix(adcpo.config.beam_angle,
    #                                 adcpo.config.beam_pattern == 'convex')

    u_inst = \
       beamvel[:,0,:] * rotmat[0, 0] +\
       beamvel[:,1,:] * rotmat[0, 1] +\
       beamvel[:,2,:] * rotmat[0, 2] +\
       beamvel[:,3,:] * rotmat[0, 3]

    v_inst = \
       beamvel[:,0,:] * rotmat[1, 0] +\
       beamvel[:,1,:] * rotmat[1, 1] +\
       beamvel[:,2,:] * rotmat[1, 2] +\
       beamvel[:,3,:] * rotmat[1, 3]

    w_inst = \
       beamvel[:,0,:] * rotmat[2, 0] +\
       beamvel[:,1,:] * rotmat[2, 1] +\
       beamvel[:,2,:] * rotmat[2, 2] +\
       beamvel[:,3,:] * rotmat[2, 3]

    errvel = \
       beamvel[:,0,:] * rotmat[3, 0] +\
       beamvel[:,1,:] * rotmat[3, 1] +\
       beamvel[:,2,:] * rotmat[3, 2] +\
       beamvel[:,3,:] * rotmat[3, 3]

    return u_inst, v_inst, w_inst, errvel

#def inst2earth(adcpo, fixed_orientation=False):
def inst2earth(u_inst, v_inst, w_inst,\
        heading_deg, pitch_deg, roll_deg,\
        orientation='down',
        heading_offset=None, declination=None,
        fixed_orientation=False):
    """Rotate velocities from the instrument to the earth frame.

    The rotation matrix is taken from the Teledyne RDI
    ADCP Coordinate Transformation manual January 2008
    """
    r = roll_deg * deg2rad
    p = np.arctan(np.tan(pitch_deg * deg2rad) * np.cos(r))
    #p = pitch_deg*deg2rad
    h = heading_deg * deg2rad

    if heading_offset is not None:
        h += heading_offset * deg2rad
    if declination is not None:
        h += declination * deg2rad
    if orientation == 'up':
        print('Adding pi to roll for up facing instrument')
        r += np.pi

    ch = np.cos(h)
    sh = np.sin(h)
    cr = np.cos(r)
    sr = np.sin(r)
    cp = np.cos(p)
    sp = np.sin(p)

        
    u = (ch * cr + sh * sp * sr) * u_inst +\
                    (sh * cp) * v_inst +\
                    (ch * sr - sh * sp * cr) * w_inst
    v = (-1 * sh * cr + ch * sp * sr) * u_inst +\
                    (ch * cp) * v_inst +\
                    (-1 * sh * sr - ch * sp * cr) * w_inst
    w = (-cp * sr) * u_inst +\
                    sp * v_inst\
                    + (cp * cr) * w_inst

    ## Multiply by the tranpose of M
    #u = (ch*cr + sh*sp*sr) *u_inst+\
    #    (ch*sp*sr - cr*sh) *v_inst+\
    #    (-cp*sr) *w_inst
    #
    #v = (cp*sh) *u_inst+\
    #    (ch*cp) *v_inst+\
    #    (sp) *w_inst

    #w = (ch*sr - cr*sh*sp) *u_inst+\
    #    (-ch*cr*sp - sh*sr) *v_inst+\
    #    (cp*cr) *w_inst

    return u, v, w

##########
# Sidelobe trimming routines
# at the moment these are set to nan, not flagged. 


##########
# Depth-bin mapping routines
def binmap(beamvel, distance,\
        pitch_deg, roll_deg,):
    """
    Map the beam data onto the appropriate depth layer
    """

    r = roll_deg * deg2rad
    p = np.arctan(np.tan(pitch_deg * deg2rad) * np.cos(r))
    #p = pitch_deg*deg2rad

    # Calculate the distance of each beam
    dist_1 = distance[..., np.newaxis] * np.cos(p)
    dist_2 = distance[..., np.newaxis] * np.cos(-p)
    dist_3 = distance[..., np.newaxis] * np.cos(r)
    dist_4 = distance[..., np.newaxis] * np.cos(-r)

    nz, nt = dist_1.shape

    beamvelnew = np.zeros_like(beamvel)
    # Build the interpolation objects
    def interp(d1, vel, dist):
        F = interp1d(d1, vel, axis=0, kind='nearest', bounds_error = False)
        return F(dist)

    print('Re-mapping beam data...')
    for ii in range(nt):
        beamvelnew[:,0,ii] = interp(dist_1[:,ii], beamvel[:,0,ii], distance)
        beamvelnew[:,1,ii] = interp(dist_2[:,ii], beamvel[:,1,ii], distance)
        beamvelnew[:,2,ii] = interp(dist_3[:,ii], beamvel[:,2,ii], distance)
        beamvelnew[:,3,ii] = interp(dist_4[:,ii], beamvel[:,3,ii], distance)

    print('Done')
    
    return beamvelnew

def _rotate_velocity_uhdas(self, beamvel,\
    heading_deg, pitch_deg, roll_deg,\
    beam_angle, beam_pattern, orientation):
    """
    Rotate from beam to compass coordinates

    Using the rotation routines from the uhdas/pycurrents library

    """

    # Create the transform objecto
    tr = Transform(angle=beam_angle, geometry=beam_pattern)

    # transform the beam to instrument coordinate
    bvel = beamvel.swapaxes(0,1).swapaxes(0,2)
    xyze = tr.xyz_to_beam(bvel)

    # rotate to east, north, up coordinates
    uvwe = rdi_xyz_enu(xyze,\
            heading_deg, pitch_deg, roll_deg,\
            orientation=orientation)

    return uvwe[:,:,0].T, uvwe[:,:,1].T, uvwe[:,:,2].T


def adcp_vars(rr):
    adcp_vars_rotation = {
        'east_vel':
            {'data':rr.ds.vel.isel(dir=0),
            'attrs':{'long_name':'Eastward water velocity',
        'coordinates':'time distance',
                'units':'m/s'},
            'dims':('distance','time')
            },
        'north_vel':
            {'data':rr.ds.vel.isel(dir=1),
            'attrs':{'long_name':'Northward water velocity',
        'coordinates':'time distance',
                'units':'m/s'},
            'dims':('distance','time')
            },
        'up_vel':
            {'data':rr.ds.vel.isel(dir=2),
            'attrs':{'long_name':'Vertical water velocity',
        'coordinates':'time distance',
                'units':'m/s'},
            'dims':('distance','time')
            },
        'errvel':
            {'data':rr.ds.vel.isel(dir=3),
            'attrs':{'long_name':'Error velocity',
        'coordinates':'time distance',
                'units':'m/s'},
            'dims':('distance','time')
            },
            }
    var_names = list(adcp_vars_rotation.keys())

    adcp_dims = {
        'distance':rr.ds['distance'].values.astype('float64'),
        'time': rr.ds['time'],
        'beam':list(range(1,5)),
    }
    
    # encoding = self.encoding
    for var in var_names:
        coords = OrderedDict()
        for dd in adcp_vars_rotation[var]['dims']:
            coords.update({dd:adcp_dims[dd]})
        V = xr.DataArray( adcp_vars_rotation[var]['data'],\
            dims=adcp_vars_rotation[var]['dims'],\
            name=var,\
            attrs = adcp_vars_rotation[var]['attrs'],\
            coords = coords            )
        rr.ds.update({var:V})

    rr.associate_qc_flag('east_vel', 'velocity')
    rr.associate_qc_flag('north_vel', 'velocity')
    rr.associate_qc_flag('up_vel', 'velocity')
    return rr