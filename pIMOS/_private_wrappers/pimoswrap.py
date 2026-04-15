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
from pIMOS.utils.file import parse_infile

# import zutils.xrwrap as xrwrap
# import zutils.xrwrap as xrwrap
import pIMOS._private_wrappers.xrwrap as xrwrap
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
    _default_attrs['pimos_nickname'] = 'pimoswrap'

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        # Ensure every subclass has its own defaults dictionary and gets a
        # deterministic nickname default based on its class name.
        if hasattr(cls, '_default_attrs') and cls._default_attrs is not None:
            cls._default_attrs = cls._default_attrs.copy()
        else:
            cls._default_attrs = cls.parent_class._default_attrs.copy()
        cls._default_attrs['pimos_nickname'] = cls.__name__

    def __init__(self):
        print('pl0 init')
        print('Initialising {}'.format(self._default_attrs['process_level']))

    @classmethod
    def blank(cls, process_level, attrs=None, verbose=False):
        """Create a wrapper instance backed by an empty in-memory Dataset.

        This is useful for filename generation or metadata workflows where no
        source file has been loaded yet.
        """

        rr = cls.__new__(cls)
        rr.verbose = verbose

        ds = xr.Dataset()
        ds.attrs = cls._default_attrs.copy()
        rr.ds = ds

        if isinstance(process_level, str):
            pl_text = process_level.strip().lower()
            if pl_text.startswith('process level'):
                process_level = pl_text.split()[-1]

        try:
            process_level_int = int(process_level)
        except Exception as error:
            raise ValueError(
                "process_level must be an integer or 'Process Level N' string."
            ) from error

        rr.update_attribute('process_level', 'Process Level {}'.format(process_level_int))

        if attrs:
            rr.update_attributes_with_dict(attrs)

        return rr

    @classmethod
    def from_netcdf(cls, infile):
        """Class-level loader for netCDF sources.

        Returns
        -------
        rr : cls
            Instantiated wrapper object.
        ds : xarray.Dataset
            Loaded dataset with load-location attrs updated.
        """

        folder, file = parse_infile(infile)
        ds = xr.open_dataset(os.path.join(folder, file))

        ds.attrs['last_load_file_name'] = file
        ds.attrs['last_load_directory'] = folder

        rr = cls(ds)
        return rr, ds
        
    @property
    def _required_attrs(self):
        """Simple function to replace the _required_attrs property of xrwrap.xrwrap to support different "Process Levels".
        
        """

        process_level = self._default_attrs['process_level']
        if process_level in ['', 0, 1]:
            _required_attrs = self.parent_class._required_attrs.copy()
            if 'pimos_nickname' not in _required_attrs:
                _required_attrs = dict.fromkeys(_required_attrs)
                _required_attrs['pimos_nickname'] = self.__class__.__name__
                _required_attrs = _required_attrs.keys()
                
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
                                                    'process_level': '',                         # This can be a default. Make it 1 when sure! 
                                                    'pimos_nickname': self.__class__.__name__,             # Default nickname follows wrapper class name.
                                                    'is_profile_data': 0}.keys()
        else:
            raise(Exception("{} is not a valid process level. You've cooked it."))

        return _required_attrs

    def assign_nickname(self, nickname):
        """
        Assign the pimos_nickname attribute. Needed for export. 
        """

        self.ds.attrs['pimos_nickname'] = nickname

#%%
