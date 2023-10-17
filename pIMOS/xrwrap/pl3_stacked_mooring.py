# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:17:42 2018

@author: 20270917
"""
#%%


import matplotlib.pyplot as plt
import cmocean as cm
import numpy as np
import xarray as xr
import matplotlib
import gsw
# import afloat.time as ztime
import pIMOS.xrwrap.pimoswrap as pimoswrap 
from iwaves.utils.fitting import fit_bmodes_linear
# from iwaves.utils.density import fit_rho, single_tanh_rho, double_tanh_rho_new


font = {'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

class_attrs = {
            'title':  'Level 3 process file made from level 2 stacked mooring',
            'source': 'pIMOS',
            'process_level': 3 
        }


def from_fv02(mooring_nc_fv02, S=35.6, **kwargs):
    """
    FUNCTION TO TAKE THE FV02 FILE AND SPIN OUT THE FV03 FILE
    """
    
    # dt_sec = kwargs.pop('dt_sec', 60)
    # start = kwargs.pop('start', None)
    # end = kwargs.pop('end', None)
    # z_method = kwargs.pop('z_method', 'z_nom')
    
    ds_fv02 = xr.open_dataset(mooring_nc_fv02)
    
    ds_fv03 = ds_fv02.copy()
    ds_fv03.attrs.update(class_attrs)

    rr = PL3_STACKED_MOORING(ds_fv03)

    return rr


class PL3_STACKED_MOORING(pimoswrap.pimoswrap):

    def __init__(self, ds, verbose=True):
        self.verbose = verbose
        if self.verbose:
            print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        self.enforce_these_attrs(class_attrs)

        # Possible functions??
        #     knockdown correct [currently written in the FV02 file but we don't call it]
        #     fill nans
        #     calc salinity
        #        1. constant
        #        2. defined T-S linear [m, b] where S = m*T + b
        #        3. use own cond data and derive T-S rel
        #     calc density
        #     set_time_slow
        #     fit modes
        #     add flags
        #     short term harmonic fit

    def knockdown_correct(self):
        """
        SHOULD BE DONE AT LEVEL 2.
        """
        raise(Exception)


    def fill_nans(self):

        print('Filling nans in Temperature')
        self.ds['temperature'] = self.ds['temperature'].interpolate_na('z_nom', fill_value="extrapolate")
        if np.sum(np.isnan(self.ds['temperature'])) > 0:
            self.ds['temperature'] = self.ds['temperature'].interpolate_na('time', fill_value="extrapolate")
        print('Done.')


    def calc_salinity(self, method='constant_sal', S=35.6):
        """
        Calculate salinity using one of a number of methods:

            1. Constant Salinity
            2. User specified T-S relationship [not implemented]
            3. Use own cond data and derive T-S rel [not implemented]

        At the moment only constant S is implemented

        Parameters
        ----------
        method: str
            Choice of "constant_sal", "user_ts", "own_cond".
        S: numeric
            Salinity to use when "constant_sal" method is selected.

        """

        print('Calculating Salinity')

        ds = self.ds

        self.ds['salinity'] = xr.zeros_like(ds.temperature) + S
        self.ds['salinity'].attrs['units'] = 'PSU'

        print('Done')


    def estimate_pressure(self):
        ds = self.ds
        site_depth = 0 # -1*ds.attrs['nominal_site_depth']
        p = site_depth - ds['z_nom'].values[...,np.newaxis]*np.ones((1, len(ds.time.values)))
        return p
    

    def calc_density(self, method='constant_sal', S=35.6):
        """
        Calculate density using GSW. Must have Temperature, Salinity fields. 
        z_nom is used for the pressure - tide and knockdown assumed small.
        Must have lat_nom and lon_nom also.
        """

        if self.verbose:
            print('Calculating Density')

        ds = self.ds

        # Determine whether to invert depth
        if ds['z_nom'][0] > ds['z_nom'][-1]:
            if self.verbose:
                print('Flipping depth var')
            ds['z_nom'] = -1 * ds['z_nom']

        SP = ds.salinity
        CT = ds.temperature

        if 'pressure' in ds.data_vars.keys():
            p  = ds.pressure # Don't worry about tides it's not sensitive enough
        else:
            p = self.estimate_pressure()
            if self.verbose:
                print('Pressure estimated from nominal depth.')
        # self.ds['pressure'] = xr.DataArray(data=p, coords=ds.coords)

        lon = ds.lon_nom.values
        lat = ds.lat_nom.values

        SA = gsw.SA_from_SP(SP, ds['z_nom'], lon, lat)

        # self.ds['sea_water_density'] = gsw.rho(SA, CT, ds['z_nom'])
        # self.ds['sea_water_density'].attrs['units'] = 'kg.m-3'

        self.ds['sea_water_potential_density'] = gsw.pot_rho_t_exact(SA, CT, ds['z_nom'], 0)
        self.ds['sea_water_potential_density'].attrs['units'] = 'kg.m-3'

        if self.verbose:
            print('Done')


    def calc_timeslow(self, slow_tstep_h=3):
        """
        Create a slow time for background quantities.
        """
        
        print('Creating slow time')

        time = self.ds.time.values

        slow_tstep = np.timedelta64(slow_tstep_h, 'h').astype('timedelta64[m]')
        
        t_start  = time[0].astype('datetime64[h]')
        t_end    = time[-1].astype('datetime64[h]')
        
        timeslow = np.arange(t_start, t_end + slow_tstep, slow_tstep)
        
        self.ds['timeslow'] = timeslow 
        # ADD ATTRIBUTES
        self.ds['timeslow'].attrs['slow_tstep_h'] = slow_tstep_h
        
        print('Done')


    def calc_bmodes(self, bar_tstep_h, density_func='double_tanh_new'):
        """
        Calculate buoyancy amplitude of internal wave modes using iwaves.fit_bmodes_linear
        """

        print('Calculating buoyancy modes')
        
        # Work out the number of suitable modes based on the instrument spacing
        # dz = signdz*2*np.diff(Z)
        ### some missing code of Raysons here
        nmodes = 4
        modes = np.arange(0,nmodes)

        # Set frequency dimensions
        # freqs = ['M2','M4','M6']
        # Nfreq = len(freqs)
        # freqvals, frqnames = harmonic_analysis.getTideFreq(freqs)

        # Set depth grid
        zmin = -280
        H = np.abs(zmin)
        zmax = 0
        Nzfit = 100
        zfit = np.linspace(zmin, zmax, Nzfit)

        rho      = self.ds['sea_water_potential_density']
        mask     = ~np.all(np.isnan(rho), axis=0)
        time     = self.ds['time']
        timeslow = self.ds['timeslow']
        Z        = self.ds['z_nom'].values # Should use z_hat if knockdown corrected

        slow_tstep_h = self.ds['timeslow'].attrs['slow_tstep_h']
        slow_tstep = np.timedelta64(slow_tstep_h, 'h').astype('timedelta64[m]')

        bar_tstep = np.timedelta64(bar_tstep_h, 'h').astype('timedelta64[m]')

        # Initiate the dataset of fitted variables
        ds_bmodes = gen_empty_bmodes_ds(time, timeslow, Nzfit+1, nmodes, zmin, zmax=0)

        
        for ii, ts in enumerate(timeslow[:-1]):

            # Get times between timesteps
            fit_ix = (time >= ts) & (time < timeslow[ii+1])

            # Get times within the background window
            t_fill = (time >= (ts + slow_tstep/2 - bar_tstep/2)) &\
                    (time < ts + slow_tstep/2 + bar_tstep/2)

            # Mask rows with any nans in them
            outx = mask & fit_ix

            # Extract background and instantaneous rho
            rhobar = np.nanmean(rho[:,t_fill], axis=1)
            rho_2fit = rho[:,outx]

            if np.sum(np.sum(outx)) > 0:
                A_t, phi, rhofit, rhofit_full, iw = \
                                fit_bmodes_linear(rho_2fit, rhobar, Z, zmin, modes,\
                                Nzfit, density_func=density_func)
            
                # Update the Dataset
                if iw.Fi.status == 0:
                    ds_bmodes['A_n'][outx,:] = A_t.T
                    ds_bmodes['rhofit'][outx,:] = rhofit_full
                    ds_bmodes['rhobar'][ii,:] = iw.rhoZ
                    ds_bmodes['N2'][ii,:] = iw.N2
                    ds_bmodes['phi'][ii,:,:] = phi.T
                    # ds_bmodes['r10'][ii,...] = r10
                    # ds_bmodes['cn'][ii,...] = c1

                    # Compute harmonics of A
        #             for ff in np.arange(0, nmodes):
        #                 An = timeseries(time[outx], A_t[ff])
        #                 amp, phs, frq, mean, yfit, yrms = An.tidefit(frqnames=freqs,\
        #                         basetime=timeslow[0])

        #                 ds_bmodes['Atide'][outx,ff] = yfit
        #                 ds_bmodes['amp'][ii,ff,:] = amp
        #                 ds_bmodes['phs'][ii,ff,:] = phs  

        self.ds['modes'] = ds_bmodes.modes
        self.ds['z_interp'] = ds_bmodes.z_interp

        for dv in list(ds_bmodes.data_vars):
            self.ds[dv] = ds_bmodes[dv]

        print('Done')

    
    def plot_background(self, warnbars=False):
        deninv_phs = check_monotonicity(self.ds, 'rhofit', 'z_interp')
        fig, ax = plt.subplots(2,1, figsize=(12,6), gridspec_kw={'hspace':0.05})
        fg = self.ds['rhobar'].T.plot(ax=ax[0], cmap=cm.cm.deep, cbar_kwargs={'pad':0.01})
        ax[0].set_xticklabels('')
        fh = self.ds['N2'].T.plot(ax=ax[1], cmap=cm.cm.solar, vmin=0., vmax=0.0004, cbar_kwargs={'pad':0.01})
        ax[1].set_xlabel('')
        ax[1].set_title('')
        fg.colorbar.set_label(label=(self.ds['rhobar'].attrs['long_name'] + '\n[' + \
                                self.ds['rhobar'].attrs['units'] + ']'), size=10)
        fh.colorbar.set_label(label=(self.ds['N2'].attrs['long_name'] + '\n[' + \
                                self.ds['N2'].attrs['units'] + ']'), size=10)
        if warnbars:
            add_warnbars(self.ds['time'], deninv_phs, ax)
        return fig, ax

        
    
def check_monotonicity(xarr, xvar, dfdim):
    den_diff = xarr[xvar].diff(dim=dfdim)
    inv_bool = np.any(den_diff < 0, axis=1)
    return inv_bool

def add_warnbars(time, warn_idx, ax):  
    # Loop thorugh warnings and axes
    for txx in time[warn_idx].values:
        for x in ax:
            # Get plot bounds
            x.axvline(txx, ymin=0, ymax=0.1, c='k', lw=0.2)
            x.axvline(txx, ymin=0.9, ymax=1., c='k', lw=0.2)

def gen_empty_bmodes_ds(time, timeslow, Nzfit, Nmode, zmin, zmax=0, freqvals=[1]):

    # Define the coordinates:
    coords = {'time':time,
                'timeslow': timeslow,
                'modes':np.arange(Nmode),
                'z_interp': np.linspace(zmin, zmax, Nzfit)[::-1], # top to bottom
                'freqs':freqvals}

    # Create the output variable as xr.DataArray objects
    A_nt = xr.DataArray(np.full((len(time), Nmode), np.nan),
            dims = ('time','modes'),
            coords = {'time':coords['time'], 'modes':coords['modes']},
            attrs = {'long_name':'Modal buoyancy amplitude',
                    'units':'m'})

    rhofit_t = xr.DataArray(np.full((len(time), Nzfit), np.nan),
            dims = ('time','z_interp'),
            coords = {'time':coords['time'], 'z_interp':coords['z_interp']},
            attrs = {'long_name':'Best-fit density',
                    'units':'kg m^-3'})

    Atide = xr.DataArray(np.full((len(time), Nmode), np.nan),
            dims = ('time','modes'),
            coords = {'time':coords['time'],'modes':coords['modes']},
            attrs = {'long_name':'Tidal fit amplitude',
                    'units':'m'})

    rho_t = xr.DataArray(np.full((len(timeslow), Nzfit), np.nan),
            dims = ('timeslow','z_interp'),
            coords = {'timeslow':coords['timeslow'], 'z_interp':coords['z_interp']},
            attrs = {'long_name':'Background density',
                    'units':'kg m^-3'})

    N2_t = xr.DataArray(np.full((len(timeslow), Nzfit), np.nan),
            dims = ('timeslow','z_interp'),
            coords = {'timeslow':coords['timeslow'], 'z_interp':coords['z_interp']},
            attrs = {'long_name':'Backgrouna buoyancy frequency squared',
                    'units':'s^-2'})

    cn_t = xr.DataArray(np.full((len(timeslow), Nmode), np.nan),
            dims = ('timeslow','modes'),
            coords = {'timeslow':coords['timeslow'], 'modes':coords['modes']},
            attrs = {'long_name':'Linear phase speed',
                    'units':'m s^-1'})

    r10_t = xr.DataArray(np.full((len(timeslow), Nmode), np.nan),
            dims = ('timeslow','modes'),
            coords = {'timeslow':coords['timeslow'], 'modes':coords['modes']},
            attrs = {'long_name':'Nonlinear steepening parameter',
                    'units':'m^-1'})

    phi_t = xr.DataArray(np.full((len(timeslow), Nmode, Nzfit), np.nan),
            dims = ('timeslow','modes','z_interp'),
            coords = {'timeslow':coords['timeslow'], 'modes':coords['modes'], 'z_interp':coords['z_interp']},
            attrs = {'long_name':'Modal structure function',
                    'units':''})

    amp_t = xr.DataArray(np.full((len(timeslow), Nmode, len(freqvals)), np.nan),
            dims = ('timeslow', 'modes', 'freqs'),
            coords = {'timeslow':coords['timeslow'],\
                    'modes':coords['modes'],\
                    'freqs':coords['freqs']},
            attrs = {'long_name':'Mode-one harmonic amplitude',
                    'units':'kg m^-3'})

    phs_t = xr.DataArray(np.full((len(timeslow), Nmode, len(freqvals)), np.nan),
            dims = ('timeslow', 'modes', 'freqs'),
            coords = {'timeslow':coords['timeslow'],\
                    'modes':coords['modes'],\
                    'freqs':coords['freqs']},
            attrs = {'long_name':'Mode-one harmonic phase',
                    'units':'radians'})

    ds = xr.Dataset({
            'A_n':A_nt,
            'rhofit':rhofit_t,
            'Atide':Atide,
            'rhobar':rho_t,
            'N2':N2_t,
            'phi':phi_t,})
#             'amp':amp_t,
#             'phs':phs_t})
            # 'r10':r10_t,
            # 'cn':cn_t})

    return ds