import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime
import numpy as np 
import pandas as pd
import glob

tind = np.arange(3000, 8000)

def stacked_scalar(ds, fh=6, kd_correct=False, tind=None, xrot=20, variable='Temperature'):
    
    plt.figure(figsize=(15, fh))
    
    if tind is None:
        tind = np.arange(len(ds.time.values))

    z = ds.z_nom
    if kd_correct:
        z = ds.z_hat[:, tind]
        
    plt.pcolor(ds.time[tind], z, ds[variable][:, tind], cmap='Spectral_r', shading='auto')
    plt.colorbar(label=variable)
    plt.ylabel('z [m]')
    plt.title(ds.attrs['site_station'])
    
    [plt.text(ds.time.values[tind][0], 
              ds.z_nom.values[n], 
              '{:.1f} m: {}'.format(ds.z_nom.values[n], ds.source.values[n]), 
              va = 'center',
              ha='left',
              fontsize=9) for n in range(0, len(ds.source.values), 2)]
    
    
    [plt.text(ds.time.values[tind][-1], 
              ds.z_nom.values[n], 
              '{:.1f} m: {}'.format(ds.z_nom.values[n], ds.source.values[n]), 
              va = 'center',
              ha='right',
              fontsize=9) for n in range(1, len(ds.source.values), 2)]

    plt.xticks(rotation=xrot)


def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))


def stacked_temp_check(ds, tempvar, dp):
    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.barh(ds['z_nom'], ds[tempvar].mean('time'),\
            height=10, color=mpl.colormaps["Spectral_r"](np.linspace(0,1,len(ds['z_nom']))))
    if not strictly_increasing(ds[tempvar].mean('time')):
        ixc = np.where(np.diff(ds[tempvar].mean('time')) < 0)[0][0]
        ax.annotate('Warning! Mean temp\nnot increasing here!',\
            (ds[tempvar].mean('time')[ixc]+0.5, ds['z_nom'][ixc]-10), fontsize=10)
    ax.set_title(dp)
    return fig
