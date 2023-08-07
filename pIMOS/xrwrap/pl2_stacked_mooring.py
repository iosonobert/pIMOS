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

# import zutils.xrwrap as xrwrap
# import zutils.stats as zstats
# import zutils.file as zfile
# import zutils.time as ztime

import afloat.time as ztime

# import pIMOS.xrwrap.pimos_pls as pls
# import pIMOS.xrwrap.seabird_37_39_56 as wrap_sbd
# import pIMOS.xrwrap.wetlabs_ntu as wetlabs_ntu
# import pIMOS.xrwrap.nortek_signature as nortek_signature
# import pIMOS.xrwrap.nortek_vector as nortek_vector 
import pIMOS.xrwrap.pimoswrap as pimoswrap 

font = {'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

class PL2_STACKED_MOORING(pimoswrap.pimoswrap):

    def __init__(self, ds):
        print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)

        class_attrs = {
            'title': 'Mooring made by stacking multipe files',
            'source': 'pIMOS' 
        }
        self.enforce_these_attrs(class_attrs)

    def knockdown_correct(rr, bottom_time, bottom_pressure, bottom_hasb):
        """
        Apply bottom correction based on another input pressure series. 

        Ignores pressure differences due to anything except BT. 

        Inputs:
            bottom_time:     Time for the bottom pressure series
            bottom_pressure: Pressure for the bottom pressure series [dB]
            bottom_hasb:     Height above bed for the bottom pressure series [m]
        """

        time_secs = ztime.seconds_since(rr.ds.time.values)
        bottom_time_secs = ztime.seconds_since(bottom_time)
        
        bottom_pressure_interp = np.interp(time_secs, bottom_time_secs, bottom_pressure)
        
        some_pressure = ~np.all(np.isnan(rr.ds.Pressure.values), axis=1)
        highest_P_z = max(rr.ds.z_nom[some_pressure].values)
        highest_P_i = np.where(rr.ds.z_nom == highest_P_z)[0][0]
        highest_P = rr.ds.Pressure[highest_P_i, :]

        pressure_diff = bottom_pressure_interp-highest_P.values
        pressure_diff_nom = highest_P_z-rr.attrs['nominal_site_depth']-bottom_hasb
        
        pressure_diff_rel = pressure_diff/pressure_diff_nom
        
        hasb = rr.ds.z_nom-rr.attrs['nominal_site_depth']
        z_hat = hasb.values[:, None]*pressure_diff_rel + rr.attrs['nominal_site_depth']
        rr.ds['z_hat'] = xr.DataArray(data=z_hat, dims=['z_nom', 'time'])

def from_fv01_archive(files, stack_variables, **kwargs):
    """
    FUNCTION TO
    NOTES: 
        - Interpolation is by nearest neighbour. A large reduction in timestep won't result in any kind of filter.
        - The first and last values of any record are set to np.nan to prevent extrap
        - Doesn't currently do any validation that the input files are in fact from the same mooring 
    
    inputs:
    
    kwargs:
        dt_sec: timestep of the stacked mooring in seconds [default=60]
        start: datatime.datetime for the start of the stacked file [default is start of first file]
        end: datatime.datetime for the end of the stacked file [default is end of first file]
        z_method: variable name used as z_coordinate in stacking [default is None in which case it tries to use file attrs]
    """
    
    dt_sec = kwargs.pop('dt_sec', 60)
    start = kwargs.pop('start', None)
    end = kwargs.pop('end', None)
    z_method = kwargs.pop('z_method', 'z_nom')
    
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
                    'nominal_instrument_height_asb',
                    'timezone',
                    ]
    attrs_that_must_be_equal = ['project',
                                'trip',
                                'site',
                                'site_station',
                                'disclaimer',
                                'nominal_latitude',
                                'nominal_longitude',
                                'nominal_site_depth',
                                'timezone']
    
    joined_attrs = {}
    for attr in attrs_to_join:
        joined_attrs[attr] = []

    stacked = {}
    
    # Used later on in validation
    first_dimension_set = None
        
    assert(type(dt_sec)==int)
    
    for i, file in enumerate(files):

        print('Loading file {} of {}'.format(i, len(files)))

        ds = xr.open_dataset(file)
        rr = pimoswrap.pimoswrap() # Just use the base class here
        rr.wrap(ds)

        for attr in attrs_to_join:
            joined_attrs[attr] += [ds.attrs[attr]]
        
        if i==0: # Just do this once. Don't need to do it at all actually
            ref_time = ds.time.values[0]
            
        print('    File starts {}'.format(ds.time.values[0]))
        time_secs = ztime.seconds_since(ds.time.values, ref_time)
        
        z_source_here = [ds.attrs['instrument_model'] + ' ' +  str(ds.attrs['instrument_serial_number'])]
        
        if i==0:
            # Initialise the time array. For now it will be clipped to the first file. 

            dt_file = ds.time.values[1] - ds.time.values[0]

            if start is None:
                start_sec = 0
            else:
                start = np.datetime64(start)
                start_sec = ztime.seconds_since(start, ref_time)
                
            if end is None:
                end_sec = time_secs[-1]
            else:
                end = np.datetime64(end)
                end_sec = ztime.seconds_since(end, ref_time)
                
            time_stacked = np.arange(start_sec, end_sec, dt_sec)

        for stack_variable in stack_variables:
            
            if stack_variable in rr.ds.data_vars.keys():
                var = rr.get_qaqc_var(stack_variable)
                var_values = var.values
                # Some validation
                if var.dims == ('distance', 'time'):
                    pass
                elif var.dims == ('time',):
                    var_values = var_values[None, :] # To make the interpolation syntax consistent for profiles and timeseries
                else:
                    print(var.dims)
                    raise(Exception('Yeahhhhhhhh nahhhhhh!!! This stacking tool has only been cobbled together for certian use cases. And unfortunately yours is not one. I recommend hacking up this code to make your own merge function'))

                # Some more validation 
                if first_dimension_set is None: # This will be None on the first pass
                    first_dimension_set = var.dims
                elif not var.dims == first_dimension_set:
                    raise(Exception("The set of variables you are trying to merge don't have trhe same dimensions. This won't work."))
                
                # KILL FIRST AND LAST VALUE
                var_values[:, 0] = np.nan
                var_values[:, -1] = np.nan
                
                # WARNING - THIS JUST DOES A NEAREST INTERP. THIS IS NOT A FILTER OR MOVING AVERAGE
                var_interp = np.array([np.interp(time_stacked, time_secs, var_row) for var_row in var_values])
                
            else:
                # Just nans if the var doesn't exist
                var_interp = np.nan*np.zeros(time_stacked.shape)
                
            if i==0:
                stacked[stack_variable] = var_interp
            else:
                stacked[stack_variable] = np.vstack([stacked[stack_variable], var_interp])


        if z_method is None:
            raise(Exception)
            z_here = ds.attrs['nominal_site_depth'] + ds.attrs['nominal_instrument_height_asb']
            z_here = np.array([z_here])[:, None]
        else:
            z_here = ds[z_method].values
            if z_here.shape==(): # If zero dimensional make one dimensional
                z_here = z_here[None] # 
            z_here = z_here[:, None] # Now make two dimensional
            
        if i==0:
            z_stacked = z_here
        else:
#             return z_stacked, z_here
            z_stacked = np.vstack([z_stacked, z_here])
            
        z_source_here = z_source_here*len(z_here)
        if i==0:
            z_source = z_source_here
        else:
            z_source += z_source_here
    
    # Make a dictionary
    stacked['time'] = np.array([ref_time + np.timedelta64(int(i), 's') for i in time_stacked])
    z_stacked = z_stacked.flatten()
    ai = np.argsort(z_stacked, axis=0)
    z_stacked = z_stacked[ai]
    stacked['z'] = z_stacked
    stacked['ai'] = ai

    for stack_variable in stack_variables:

        stacked[stack_variable] = stacked[stack_variable][ai, :]
        
    #########################
    ## Make the xr DataSet ##
    #########################
    ds_stacked = xr.Dataset(attrs={})
    encoding = {}
    coords = OrderedDict()
    coords.update({z_method: z_stacked})
    coords.update({'time': stacked['time']})
    
    for stack_variable in stack_variables:
        
        V = xr.DataArray( stacked[stack_variable],\
            dims=coords,\
            name=stack_variable,\
            attrs = {},\
            coords = coords
        )

        ds_stacked.update({stack_variable:V})
        
        encoding.update({stack_variable:{'zlib':True,'_FillValue':-999999.}})
        
    z_source = np.array(z_source)[ai]
    V = xr.DataArray( z_source,\
            dims={z_method: z_stacked},\
            name='source',\
            attrs = {},\
            coords = {z_method: z_stacked}
        )
    ds_stacked.update({'source':V})
    
    rr = PL2_STACKED_MOORING(ds_stacked)

    ####################
    ## ATTRIBUTE WORK ##
    ####################
    # FIRST THE ATTRIBUTES WHICH MUST BE EQUAL
    for attr in attrs_that_must_be_equal:
        # split_attr = joined_attrs[attr].split(';')              # Split back to list
        split_attr = joined_attrs[attr]          
        if not np.all(np.array(split_attr) == split_attr[0]):
            if isinstance(split_attr[0], float) &\
                (np.all(np.around(split_attr,5) == np.around(split_attr[0],5))):
                print('{} is a little off, but close enough'.format(attr))
            # else:
            #     print(split_attr)
            #     raise(Exception('All files must have the same {}'.format(attr)))
        rr.ds.attrs[attr] = split_attr[0]                       # Assign just the first value
        joined_attrs.pop(attr)                                  # Remove from dict

    # NEXT THE ATTRIBUTES WHICH CAN BE JOINED
    rr.ds.attrs['source'] = ';'.join(files)
    for attr in joined_attrs.keys():
        # print(joined_attrs[attr])
        joined_attrs[attr] = [str(i) for i in joined_attrs[attr]]
        rr.ds.attrs[attr] = ';'.join(joined_attrs[attr])

    return rr


# Function to find all files in archive folder for a list of selected instruments
def get_files(archive_dir, deployment, mooring, instruments):
    files = []
    for instrument in instruments:
        f_format = '*[[]{}[]]*[[]{}[]]*[[]{}[]]*.nc'.format(deployment, mooring, instrument)
        files += glob.glob(os.path.join(archive_dir, f_format))
    
    if len(files) == 0:
        raise(Exception("No files found"))
        
    return files


def check_files(files, var):
    files_in = np.full(len(files), False)
    for xx, file in enumerate(files):
        ds = xr.open_dataset(file)
        if var in ds.data_vars.keys():
            files_in[xx] = True
    return files_in
    
