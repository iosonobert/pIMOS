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

import zutils.xrwrap as xrwrap
import zutils.stats as zstats
import zutils.file as zfile

font = {'weight' : 'normal',
        'size'   : 15}
matplotlib.rc('font', **font)

class PIMOS_PL0(xrwrap.xrwrap):

    parent_class = xrwrap.xrwrap
    _default_attrs = parent_class._default_attrs.copy()
    _default_attrs['process_level'] = 'Process Level 0'

    def __init__(self):
        print('pl0 init')
        print('Initialising {}'.format(self._default_attrs['process_level']))
        

class PIMOS_PL1(PIMOS_PL0):

    parent_class = PIMOS_PL0
    _default_attrs = parent_class._default_attrs.copy()
    _default_attrs['process_level'] = 'Process Level 1'
    
    def __init__(self):
        print('pl1 init')
        self.parent_class.__init__(self)
        
        pass
        
class PIMOS_PL2(PIMOS_PL1):
        
    parent_class = PIMOS_PL1
    _default_attrs = parent_class._default_attrs.copy()
    _default_attrs['process_level'] = 'Process Level 2'

    _required_attrs = {
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

    def __init__(self):
        print('pl2 init')

        self.parent_class.__init__(self)
        
class PIMOS_PL3(PIMOS_PL2):

    parent_class = PIMOS_PL2
    _default_attrs = parent_class._default_attrs.copy()
    _default_attrs['process_level'] = 'Process Level 3'

    def __init__(self):

        self.parent_class.__init__(self)

        pass

class PIMOS_PL4(PIMOS_PL3):

    parent_class = PIMOS_PL3
    _default_attrs = parent_class._default_attrs.copy()
    _default_attrs['process_level'] = 'Process Level 4'

    def __init__(self):

        self.parent_class.__init__(self)

        pass
        
#%%
