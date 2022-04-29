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
import glob

import importlib 
from collections import OrderedDict

import zutils.xrwrap as xrwrap
import zutils.stats as zstats
import zutils.file as zfile
import zutils.time as ztime

import pIMOS.xrwrap.pimos_pls as pls
import pIMOS.xrwrap.seabird_37_39_56 as wrap_sbd
import pIMOS.xrwrap.wetlabs_ntu as wetlabs_ntu
import pIMOS.xrwrap.nortek_signature as nortek_signature
import pIMOS.xrwrap.nortek_vector as nortek_vector 
import pIMOS.xrwrap.pimoswrap as pimoswrap 

font = {'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

class PL2_STACKED_TRANSECT(pimoswrap.pimoswrap):

    def __init__(self, ds):
        print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        class_attrs = {
            'title': 'Transect made by stacking multipe files',
            'source': 'pIMOS' 
        }
        self.enforce_these_attrs(class_attrs)


def from_fv01_archive_slow(files, stack_variables, znew, **kwargs):
    """
    Very slow as it uses a slow interp code and it interpolates every variable twice. 
    
    Set include_funky to False for faster interp
    """
    
    dt_sec = kwargs.pop('dt_sec', 60)
    start = kwargs.pop('start', None)
    end = kwargs.pop('end', None)
    z_method = kwargs.pop('z_method', 'z_nom')

    include_funky = kwargs.pop('include_funky', True)

    attrs_to_join = ['project',
                    'trip',
                    'trip_deployed',
                    'site',
                    'site_station',
                    'instrument_make',
                    'instrument_model',
                    'instrument_serial_number',
                    'raw_file_name',
                    'raw_file_directory',
                    'raw_file_attributes',
                    'disclaimer',
                    'nominal_latitude',
                    'nominal_longitude',
                    'nominal_site_depth',
                    'timezone',
                    'is_profile_data'
                    ]
    attrs_that_must_be_equal = ['project',
                                'trip',
                                'site',
                                'instrument_make',
                                'instrument_model',
                                'instrument_serial_number',
                                'site_station',
                                'disclaimer',
                                'timezone',
                                'is_profile_data']
    
    joined_attrs = {}
    for attr in attrs_to_join:
        joined_attrs[attr] = []

    output = {}

    nz = len(znew)
    
    # Used later on in validation
    first_dimension_set = None
        

    for i, file in enumerate(files):

        print('Loading file {} of {}'.format(i, len(files)))
        
        ds = xr.open_dataset(file)
        rr = xrwrap.xrwrap() # Just use the base class here
        rr.wrap(ds)

        for attr in attrs_to_join:
            joined_attrs[attr] += [ds.attrs[attr]]

        z = ds.Depth.values
        
        # Funky interp scheme 
        znewi = np.linspace(znew[0], max(z), nz) 
        print(max(z))
        
        for stack_variable in stack_variables:

            data = ds[stack_variable].values
            datanew = slow_interp(znew, z, data, l=2)

            if i == 0:
                output[stack_variable] = datanew[:, None]

            else:
                output[stack_variable] = np.hstack((output[stack_variable], datanew[:, None]))

            if include_funky:
                # Funky interp scheme 
                datanewi = slow_interp(znewi, z, data, l=2)
                if i == 0:
                    output[stack_variable+'_i'] = datanewi[:, None]
                else:
                    output[stack_variable+'_i'] = np.hstack((output[stack_variable+'_i'], datanewi[:, None]))
                
        latitude = np.mean(rr.ds.latitude.values)
        longitude = np.mean(rr.ds.longitude.values)
        time = ztime.datetime64_mean(rr.ds.datenum.values)
        
        if i == 0:
            latitudes = np.array([latitude])[:, None]
            longitudes = np.array([longitude])[:, None]
            times = time[:, None]
            znewis = znewi[:, None]
        else:
            latitudes = np.hstack((latitudes, np.array([latitude])[:, None]))
            longitudes = np.hstack((longitudes, np.array([longitude])[:, None]))
            times = np.hstack((times, time[:, None]))
            znewis = np.hstack((znewis, znewi[:, None]))

    output['latitude'] = latitudes
    output['longitude'] = longitudes
    output['time'] = times
    output['z'] = znew
    
    if include_funky:
        output['z_i'] = znewis
    
    # Sort
    ai = np.argsort(output['time'])
    ai = ai .flatten()
    output['latitude'] = output['latitude'][0, ai]
    output['longitude'] = output['longitude'][0, ai]
    output['time'] = output['time'][0, ai]
    
    if include_funky:
        output['z_i'] = output['z_i'][:, ai]

    for stack_variable in stack_variables:

        output[stack_variable] = output[stack_variable][:, ai]
        if include_funky:
            output[stack_variable+'_i'] = output[stack_variable+'_i'][:, ai]
        
    # Add a distance variable
    dlon = output['longitude'] - output['longitude'][0]
    dlat = output['latitude'] - output['latitude'][0]
    dist_lon = dlon*111310
    dist_lat = dlat*111310*np.cos(np.mean(output['latitude']*np.pi/180)) # reference

    output['dist_approx'] = np.sqrt(dist_lat**2 + dist_lon**2)

    #########################
    ## Make the xr DataSet ##
    #########################
    ds_stacked = xr.Dataset(attrs={})
    encoding = {}
    coords = OrderedDict()
    coords.update({'z': output['z']})
    coords.update({'time': output['time']})

    if include_funky:
        coordsi = OrderedDict()
        coordsi.update({'z_i': output['z_i']})
        coordsi.update({'time': output['time']})

        print(output['z_i'].shape)
        V = xr.DataArray( output['z_i'],\
            dims=coords,\
            name='z_i',\
            attrs = {},\
            coords = coords
        )
        ds_stacked.update({'z_i': V})

    V = xr.DataArray( output['dist_approx'],\
            dims={'time': output['time']},\
            name='dist_approx',\
            attrs = {},\
            coords = {'time': output['time']}
        )
    ds_stacked.update({'dist_approx': V})
    V = xr.DataArray( output['longitude'],\
            dims={'time': output['time']},\
            name='longitude',\
            attrs = {},\
            coords = {'time': output['time']}
        )
    ds_stacked.update({'longitude': V})
    V = xr.DataArray( output['latitude'],\
            dims={'time': output['time']},\
            name='latitude',\
            attrs = {},\
            coords = {'time': output['time']}
        )
    ds_stacked.update({'latitude': V})

    
    for stack_variable in stack_variables:

        V = xr.DataArray( output[stack_variable],\
            dims=coords,\
            name=stack_variable,\
            attrs = {},\
            coords = coords
        )

        ds_stacked.update({stack_variable:V})

        if include_funky:

            V = xr.DataArray( output[stack_variable+'_i'],\
            dims=coords,\
            name=stack_variable+'_i',\
            attrs = {},\
            coords = coords
        )

        ds_stacked.update({stack_variable+'_i':V})

        encoding.update({stack_variable:{'zlib':True,'_FillValue':-999999.}})

    # z_source = np.array(z_source)[ai]
    # V = xr.DataArray( z_source,\
    #         dims={z_method: z_stacked},\
    #         name='source',\
    #         attrs = {},\
    #         coords = {z_method: z_stacked}
    #     )
    # ds_stacked.update({'source':V})

    rr = PL2_STACKED_TRANSECT(ds_stacked)

    ####################
    ## ATTRIBUTE WORK ##
    ####################
    # FIRST THE ATTRIBUTES WHICH MUST BE EQUAL
    for attr in attrs_that_must_be_equal:
        # split_attr = joined_attrs[attr].split(';')              # Split back to list
        split_attr = joined_attrs[attr]
        if not np.all(np.array(split_attr) == split_attr[0]):
            print(split_attr)
            raise(Exception('All files must have the same {}'.format(attr)))
        rr.ds.attrs[attr] = split_attr[0]                       # Assign just the first value
        joined_attrs.pop(attr)                                  # Remove from dict

    # NEXT THE ATTRIBUTES WHICH CAN BE JOINED
    rr.ds.attrs['source'] = ';'.join(files)
    for attr in joined_attrs.keys():
        # print(joined_attrs[attr])
        joined_attrs[attr] = [str(i) for i in joined_attrs[attr]]
        rr.ds.attrs[attr] = ';'.join(joined_attrs[attr])

    return rr

def slow_interp(znew, z, data, l=2):
    """
    Slow loopy interp function I'm gonna use.
    
    """
    def get_weights(znew_, z, l):
        d = (znew_ - z)
        w = np.exp(-(d**2)/(l**2))

        return w

    datanew = np.ones((len(znew),))
    for i in np.arange(len(znew)):

        w = get_weights(znew[i], z, l=l)
        datanew[i] = sum(data*w)/sum(w)
          
    return datanew
