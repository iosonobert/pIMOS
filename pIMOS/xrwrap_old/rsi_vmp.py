import sys, os, scipy.io
import xarray as xr
import zutils.xrwrap as xrwrap
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import datetime

###############################################
# These dependencies could be loader specific #
###############################################
import pyODAS.read_p

def from_netcdf(infile):
    """
    Pass straight to the main xrwrap load method.
    """

    folder, file = xrwrap.parse_infile(infile)

    ds = xr.open_dataset(os.path.join(folder, file))

    ds.attrs['last_load_file_name']      = file
    ds.attrs['last_load_directory']      = folder

    print(ds)
    
    rr = RSI_VMP(ds)

    return rr, ds

def from_pfile(pfile, ini_file):
    """
    Reads RSI *.p file and return xarray dataset with the xrray accessor extension. 
    """

    pdict = pyODAS.read_p.main(pfile, ini_file=ini_file)
    
    ds = xr.Dataset()
    
    dict_copy = pdict.copy()

    ds.attrs['raw_file_name']      = os.path.split(pfile)[1]
    ds.attrs['raw_file_directory'] = os.path.split(pfile)[0]

    ds.attrs['Year']         = dict_copy.pop('year')
    ds.attrs['Month']        = dict_copy.pop('month')
    ds.attrs['Day']          = dict_copy.pop('day')
    ds.attrs['Hour']         = dict_copy.pop('hour')
    ds.attrs['Minute']       = dict_copy.pop('minute')
    ds.attrs['Second']       = dict_copy.pop('second')
    ds.attrs['Millisecond']  = dict_copy.pop('millisecond')
    
    ds.attrs['datetime']     = datetime.datetime(
                                ds.attrs['Year'],
                                ds.attrs['Month'],
                                ds.attrs['Day'],
                                ds.attrs['Hour'],
                                ds.attrs['Minute'],
                                ds.attrs['Second'],
                                ds.attrs['Millisecond'],
                                    )
    
    ds.attrs['config_string'] = dict_copy.pop('config_string')

    ds.attrs['fs_fast'] = dict_copy.pop('fs_fast')
    ds.attrs['fs_slow'] = dict_copy.pop('fs_slow')

    t_slow = dict_copy.pop('t_slow').flatten()
    t_fast = dict_copy.pop('t_fast').flatten()
    
    t_slow = [ds.attrs['datetime'] + datetime.timedelta(seconds=x) for x in t_slow]
    t_fast = [ds.attrs['datetime'] + datetime.timedelta(seconds=x) for x in t_fast]

    ds = ds.assign_coords(time_slow = t_slow)
    ds = ds.assign_coords(time_fast = t_fast)

    dict_copy.pop('header') # Don't really know what to do with this just yet
    dict_copy.pop('channels') # Don't really know what to do with this just yet
    dict_copy.pop('all_indices') # Don't really know what to do with this just yet
    dict_copy.pop('all_vectors') # Don't really know what to do with this just yet
    dict_copy.pop('all_times') # Don't really know what to do with this just yet

    dict_copy.pop('Gnd') # Don't really know what to do with this just yet

    # The rest should be data variables
    keys = dict_copy.keys()
    for key in keys:

        if 'data_physical' in dict_copy[key].keys():
            setname = 'data_physical'
        else:
            setname = 'data'
        if len(dict_copy[key][setname])==len(ds.time_fast):
            ds[key] = xr.DataArray(dict_copy[key][setname].flatten(), dims=["time_fast"])
        elif len(dict_copy[key][setname])==len(ds.time_slow):
            ds[key] = xr.DataArray(dict_copy[key][setname].flatten(), dims=["time_slow"])
        else:
            print('Not sure what this is: {}!'.format(key))
            
    ###############################################
    # Initialise the accessor
    #   Needs only be done once
    ###############################################
    rr = RSI_VMP(ds)
        
    return rr, ds

##########################
# Actual xarray wrap #
##########################
class RSI_VMP(xrwrap.xrwrap):
    """
    Accessor for the RSI VMP. 
    """
    def __init__(self, ds):
        
        print('Initialising accessor.')
        # self._obj = ds
        # self.wrap(ds) # XRWRAP compatibility
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)
        self.update_global_attrs()

        pass
    
    def update_global_attrs(self):
        """
        Each wrapper should overload this function
        """

        # CF Compliance
        print('Updating attributes function of the class.')

        class_attrs = {
            'title': 'Measured data from a RSI VMP read from .P file',
            'source': 'RSI VMP'
        }
        
        print(self.attrs)
        # source = 
        for attr in class_attrs.keys():
            if self.attrs[attr] == '':
                self.update_attribute(attr, class_attrs[attr])
            else:
                assert(self.attrs[attr]==class_attrs[attr])        