import os
import xarray as xr
import matplotlib.pyplot as plt
import datetime
import numpy as np 
import pandas as pd
import glob

tind = np.arange(3000, 8000)

def stacked_scalar(ds, fh=6, kd_correct=False, tind=None, xrot=20):
    plt.figure(figsize=(15, fh))
    
    if tind is None:
        tind = np.arange(len(ds.time.values))

    z = ds.z_nom
    if kd_correct:
        z = ds.z_hat[:, tind]
        
    plt.pcolor(ds.time[tind], z, ds.Temperature[:, tind], cmap='Spectral_r', shading='auto')
    plt.colorbar(label='Temperature')
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
