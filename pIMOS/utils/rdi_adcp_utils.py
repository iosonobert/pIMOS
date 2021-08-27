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

            dsnew.update({'u':u, 'v':v, 'w':w,'errvel':evel})
            
    else:
        dsnew.update({'v' : resample( ds['v'], dtavg, axis=1) })
        dsnew.update({'u' : resample( ds['u'], dtavg, axis=1) })
        dsnew.update({'w' : resample( ds['w'], dtavg, axis=1) })

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

    ds = self.ds

    evel = ds['errvel']

    mask = evel.values > thresh

    ds['u'].values[mask] = np.nan
    ds['v'].values[mask] = np.nan
    ds['w'].values[mask] = np.nan

    self.ds.attrs.update({'QAQC Error Velocity Maximum Threshold':thresh})


def qaqc_tilt(self, cutoff_angle):
    """
    Mask u/v/w arrays when the instrument tilt goes too far
    """
    ds = self.ds

    tilt = self._calc_tilt(ds['pitch'].values,\
            ds['roll'].values)

    flag_angle = tilt > cutoff_angle

    # Create a 2d mask array
    nz = ds.distance.shape
    mask = flag_angle[np.newaxis,:].repeat(nz,0)

    ds['u'].values[mask] = np.nan
    ds['v'].values[mask] = np.nan
    ds['w'].values[mask] = np.nan

    self.ds.attrs.update({'QAQC Maximum Tilt Angle':cutoff_angle})

def qaqc_depth(self, max_depth, orientation=None, P=None, variables=['u', 'v', 'w']):
    """
    Mask out regions outside of the maximum depth
    """

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

def sidelobe_trim(self, trim=1.5, P='depth', variables=['u', 'v', 'w'], trim_units='native'):
    """
    This actually sets these variables to nan. It does not flag them. 
    """

    sidelobe_trim_calc(self, trim=trim, P=P, trim_units=trim_units)
    for variable in variables:
        self.ds[variable].values[self.ds['sidelobe_blank'].values] = np.nan

def sidelobe_trim_calc(self, trim=1.5, P='depth', trim_units='native'):
    """
    inputs:
        P - the variable to trim by. Can be numeric or the name of a variable [e.g. pressure] with the appropriate dimensions.
               -- for upward facing instruments the 'pressure' variable is a good choice, or a nominal pressure if the pressure sensor failed.  
               -- for downward facing instruments the instrument height above the bed would be a good choice, or a nominal height if this is unknown.
        """

    # Make the blanking variables
    if not 'sidelobe_last_good' in self.ds.data_vars.keys():
        self.ds['sidelobe_last_good'] = xr.DataArray(
            data=np.zeros_like(self.ds.mpltime),
            dims=self.ds.mpltime.dims,
            coords=self.ds.mpltime.coords,
            attrs=dict())
        
    if not 'sidelobe_blank' in self.ds.data_vars.keys():
        self.ds['sidelobe_blank'] = xr.DataArray(
            data=np.zeros_like(self.ds.u),
            dims=self.ds.u.dims,
            coords=self.ds.u.coords,
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

    for i in np.arange(0, 10):
        print('USE NP.MATMUL')

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
        r += np.pi

    ch = np.cos(h)
    sh = np.sin(h)
    cr = np.cos(r)
    sr = np.sin(r)
    cp = np.cos(p)
    sp = np.sin(p)

    for i in np.arange(0, 10):
        print('USE NP.MATMUL')
        
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
