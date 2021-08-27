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
    
    raise(NotImplementedError)

#######
# QA/QC methods
#######
def qaqc_pgood(self, thresh):
    """
    Mask u/v/w based on percentage good information

    Uses minimum of all four beams
    """

    raise(NotImplementedError)

def qaqc_corr(self, thresh):
    """
    Mask u/v/w based on correlation

    Ensure that the minimum of three beams is above threshold

    Use beams 1-3 for now as these are used for the transformation

    TODO: Use the correlation to indicate which beams to use for transformation
    """

    raise(NotImplementedError)

    # ds = self.ds
    # corr = ds['corr'][:,0:3,:].min(axis=1)

    # mask = corr.values < thresh

    # ds['u'].values[mask] = np.nan
    # ds['v'].values[mask] = np.nan
    # ds['w'].values[mask] = np.nan

    # self.ds.attrs.update({'QAQC Correlation Minimum Threshold':thresh})

def qaqc_echo(self, thresh):
    """
    Mask u/v/w based on echo intensity

    Ensure that all four beams are above threshold
    """

    raise(NotImplementedError)

    # ds = self.ds
    # echo = ds['echo'][:,:,:].min(axis=1)

    # mask = echo.values < thresh

    # ds['u'].values[mask] = np.nan
    # ds['v'].values[mask] = np.nan
    # ds['w'].values[mask] = np.nan

    # self.ds.attrs.update({'QAQC Echo Intensity Minimum Threshold':thresh})

def qaqc_errvel(self, thresh):
    """
    Mask u/v/w based on error velocity

    """

    raise(NotImplementedError)

    # ds = self.ds

    # evel = ds['errvel']

    # mask = evel.values > thresh

    # ds['u'].values[mask] = np.nan
    # ds['v'].values[mask] = np.nan
    # ds['w'].values[mask] = np.nan

    # self.ds.attrs.update({'QAQC Error Velocity Maximum Threshold':thresh})

def qaqc_tilt(self, cutoff_angle):
    """
    Mask u/v/w arrays when the instrument tilt goes too far
    """

    raise(NotImplementedError)
    
    # ds = self.ds

    # tilt = self._calc_tilt(ds['pitch'].values,\
    #         ds['roll'].values)

    # flag_angle = tilt > cutoff_angle

    # # Create a 2d mask array
    # nz = ds.distance.shape
    # mask = flag_angle[np.newaxis,:].repeat(nz,0)

    # ds['u'].values[mask] = np.nan
    # ds['v'].values[mask] = np.nan
    # ds['w'].values[mask] = np.nan

    # self.ds.attrs.update({'QAQC Maximum Tilt Angle':cutoff_angle})

def qaqc_depth(self, max_depth, orientation=None, P=None, variables=['u', 'v', 'w']):
    """
    Mask out regions outside of the maximum depth
    """

    raise(NotImplementedError)

def sidelobe_trim(self, trim=1.5, P='depth', variables=['u', 'v', 'w'], trim_units='native'):

    raise(NotImplementedError("Take this from RDI ADCP and make generic"))

#########
# Instrument rotation routines
#########
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
    
    raise(NotImplementedError)


"""
Nortek notes on transforms
There are a couple of Matlab codes in there too.  

https://support.nortekgroup.com/hc/en-us/articles/360029820971-How-is-a-coordinate-transformation-done-

https://nortek.zendesk.com/attachments/token/g3Xal028bJkYclRph8hRdlIR2/?name=signatureAD2CP_beam2xyz_enu.m
"""
def beam2inst(beamvel, T):
    """Rotate velocities from beam to instrument coordinates.

    4 beam inst to xyzz
    
    """
    nc = beamvel.shape[1]
    XYZZ = np.zeros_like(beamvel)
    
    for i in np.arange(0, nc):
        XYZZ[:, i, :] = np.matmul(T, beamvel[0:4, i, :])
        # XYZ = np.matmul(T, T)

    four_beam_error_velocity = XYZZ[2, :, :]-XYZZ[3, :, :]
    XYZZ[2, :, :] = (XYZZ[2, :, :]+XYZZ[3, :, :])/2

    return XYZZ

def inst2earth(XYZ, heading_deg, pitch_deg, roll_deg,\
        orientation='down',
        heading_offset=None, declination=None,
        fixed_orientation=False):

    
    if orientation.lower()=='down':
        raise(NotImplementedError)
    if orientation.lower()=='up':
        pass
    
    nc = XYZ.shape[1]
    
    heading_deg = heading_deg-90
    hh = heading_deg * deg2rad
    pp = pitch_deg * deg2rad
    rr = roll_deg * deg2rad
    
    if fixed_orientation: # Speed up for fixed instrument
        hh = np.mean(hh)
        pp = np.mean(pp)
        rr = np.mean(rr)
    else:
        raise(NotImplementedError)

    # % Make heading matrix
    H = np.array([[np.cos(hh), np.sin(hh), 0],[ -np.sin(hh), np.cos(hh), 0],[0, 0, 1]])

    # % Make tilt matrix
    P = np.array([[np.cos(pp), -np.sin(pp)*np.sin(rr), -np.cos(rr)*np.sin(pp)],
         [0, np.cos(rr), -np.sin(rr)],  
         [np.sin(pp), np.sin(rr)*np.cos(pp),  np.cos(pp)*np.cos(rr)]])

    # % Make resulting transformation matrix
    xyz2enu = np.matmul(H, P)

    xyz2enu.shape
    
    ENU = np.zeros_like(XYZ[:, :, :])

    for i in np.arange(0, nc):
        ENU[0:3, i, :] = np.matmul(xyz2enu, XYZ[:, i, :])
        # XYZ = np.matmul(T, T)
    
    return ENU


##########
# Depth-bin mapping routines
def binmap(beamvel, distance,\
        pitch_deg, roll_deg,):
    """
    Map the beam data onto the appropriate depth layer
    """

    raise(NotImplementedError)

class NotImplementedError(Exception):
    pass