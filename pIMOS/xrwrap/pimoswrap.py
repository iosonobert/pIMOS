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

import importlib 

# import zutils.xrwrap as xrwrap
# import zutils.xrwrap as xrwrap
import afloat.utils.xrwrap as xrwrap
# import zutils.stats as zstats
# import zutils.file as zfile

font = {'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

class pimoswrap(xrwrap.xrwrap):
    """
    Subclassing the xrwrap.xrwrap class to give some added pIMOS specific functionality. 
    """

    parent_class = xrwrap.xrwrap
    _default_attrs = parent_class._default_attrs.copy()
    _default_attrs['process_level'] = 'Process Level 0'

    def __init__(self):
        print('pl0 init')
        print('Initialising {}'.format(self._default_attrs['process_level']))
        
    @property
    def _required_attrs(self):
        """Simple function to replace the _required_attrs property of xrwrap.xrwrap to support different "Process Levels".
        
        """

        process_level = self._default_attrs['process_level']
        if process_level in ['', 0, 1]:
            _required_attrs = parent_class._required_attrs.copy()
        elif process_level in [2, 3, 4]:
            _required_attrs = _required_attrs = {
                                                    'title': '', 
                                                    'institution': 'The University of Western Australia', 
                                                    'institution_division': 'Ocean Dynamics', 
                                                    'source': '', 
                                                    'project': '', 
                                                    'history': '', 
                                                    'references': '', 
                                                    'comment': '', 
                                                    'Conventions': 'CF-1.7', 
                                                    'site': '', 
                                                    'site_station': '', 
                                                    'last_export_file_name': '',                # When a netcdf is loaded, this should be cleared
                                                    'last_export_directory': '',                # When a netcdf is loaded, this should be cleared
                                                    'last_load_file_name': '',                  # When a netcdf is loaded, this should be overwritten with the load_name
                                                    'last_load_directory': '',                  # When a netcdf is loaded, this should be overwritten with the load_name
                                                    'outfile_append': '', 
                                                    'disclaimer': '',
                                                    'nominal_latitude': '',
                                                    'nominal_longitude': '',
                                                    'nominal_site_depth': '',
                                                    'timezone': '',
                                                    'process_level': '',
                                                    'is_profile_data': 0}.keys()
        else:
            raise(Exception("{} is not a valid process level. You've cooked it."))

        return _required_attrs

#%%
