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
from matplotlib.dates import num2date, date2num
import matplotlib
import datetime
import os 
import pdb 

from pIMOS.utils import qc_conventions
from pIMOS.utils import file as zfile

default_attrs = {
    'title': '', 
    'institution': 'The University of Western Australia', 
    'institution_division': 'Ocean Dynamics', 
    'source': '', 
    'project': '', 
    'history': '', 
    'references': '', 
    'comment': '', 
    'Conventions': 'CF-1.7', 
    'trip': '', 
    'trip_deployed': '', 
    'site': '', 
    'site_station': '', 
    'instrument_make': '',  
    'instrument_model': '',
    'instrument_serial_number': '',
    'raw_file_name': '',
    'raw_file_directory': '',
    'raw_file_attributes': '',                  # This is a string which holds the attributes of any raw file taken in by the. 
    'last_export_file_name': '',                # When a netcdf is loaded, this should be cleared
    'last_export_directory': '',                # When a netcdf is loaded, this should be cleared
    'last_load_file_name': '',                  # When a netcdf is loaded, this should be overwritten with the load_name
    'last_load_directory': '',                  # When a netcdf is loaded, this should be overwritten with the load_name
    'outfile_append': '', 
    'disclaimer': '',
    'nominal_latitude': '',
    'nominal_longitude': '',
    'nominal_site_depth': '',
    'pressure_sensor_height_asb': '',
    'nominal_instrument_height_asb': '',
    'nominal_instrument_orientation': '',
    'timezone': '',
    'process_level': '',
    'is_profile_data': 0} 

def parse_infile(infile, verbose=True):
    """
    Take a string ( folder/file ), or 2 element list  ( [folder, file] ), and return the folder and file.
    """

    if type(infile)==str:
        if verbose:
            print('Infile is a string')
        folder, file = os.path.split(infile)
    elif type(infile) in [list, tuple]:
        if not len(infile) == 2:
            raise(Exception('The infile must be a string or a length 2 sequence'))
        else:
            folder, file = infile
    else:
        raise(Exception('The infile must be a string or a length 2 sequence'))
        
    return folder, file
    
def _from_netcdf(infile, classhandler):
    """
    Pass straight to the main xrwrap load method.

    Inputs:
        - infile
        - classhandler - the xrwrap class to open the file
    """
   
    folder, file = parse_infile(infile)

    ds = xr.open_dataset(os.path.join(folder, file))

    ds.attrs['last_load_file_name']      = file
    ds.attrs['last_load_directory']      = folder

    # print(ds)
    
    rr = classhandler(ds)

    return rr, ds

class xrwrap():
    
    ######
    # SHOULD MINIMISE THE NUMBER OF PROPERTIES HERE, BECAUSE THESE DON'T EXPORT TO NETCDF OR LOAD FROM NETCDF 
    ######
    # folder = ''
    # file_ = ''

    # Certain irreversible operations will lead to modification of the output file name to prevent overwrite. 
    auto_file_append = '' 
    
    so = None
    eo = None

    first_good = None
    last_good = None

    _qc_conv = qc_conventions.base()
    _default_attrs = default_attrs
    _required_attrs = default_attrs.keys()
    _outname = None
    
    # attrs = default_attrs

    default_user = 'UWA'

    verbose = False

    @property
    def _obj(self):
        """
        Little trick that may or may not work for xarray accessors.
        """
        return self.ds

    @property
    def fullpath(self, i=0):
        
        raise(Exception('Don''t use this fullpath property anymore. Switching to a more detailed system involving dataset attributes. '))

        if type(self.file_) == list:
            return '{folder}/{file_}'.format(folder=self.folder, file_=self.file_[i])
        else:
            return '{folder}/{file_}'.format(folder=self.folder, file_=self.file_)

    @property
    def attrs(self):
        return self.ds.attrs

    def generate_outname(self, keep_name='False', fileappend=''):
        """
        Function to autogenerate a file name. Should not do this continuously.
        """
        use_profiler_covention = False
        if self.is_profile_data or use_profiler_covention:
            convention = ['project',
                            'trip',
                            'site_station',
                            'instrument_model',
                            'instrument_serial_number',
                            'raw_file_name',
                            'outfile_append']
        else:

            convention = ['project',
                            'trip',
                            'site_station',
                            'nominal_instrument_height_asb',
                            'instrument_model',
                            'instrument_serial_number',
                            'outfile_append']

        outname = naming_conv(self.attrs, convention=convention)
        return outname

    # @property
    def outname(self, fileappend=''):
        """
        generate an outname only if outname is empty.
        """
        if self._outname is None:
            self._outname = self.generate_outname(fileappend='')

        return self._outname

    def load(self, infile):
        """
        Initialise from netcdf rather than from a raw input file. 
        """

        raise(Exception("Don't load like this anymore. Load from outside the actual object."))

        folder, file = parse_infile(infile)

        nc_file = os.path.join(folder, file)
        
        ds = xr.open_dataset(nc_file)

        ds.attrs['last_load_file_name']      = file
        ds.attrs['last_load_directory']      = folder

        self.store_raw_file_attributes(ds)
        self.ds = ds

        # self.attrs = ds.attrs
        # self.ds.attrs = {}

        print('Loaded from NC. Class attributes taken from NC. NC attrs cleared.')

        pass

    def wrap(self, ds):
        """
        Initialise from an xarray dataset. Got to be carful here that the dataset is of the exact format or lots of errors will be thrown. Could used this if 2 NC files
        are merged outside of here and then brought back in. That's the only usage I can immagine.  
        """

        self.ds = ds
        if False: # Still not sure the best point to do this.
            self.store_raw_file_attributes(ds)

        print('Wrapped an existing xarray dataset. Class attributes taken from the dataset. Dataset attrs cleared.')
        
    def export(self, naming_method=None, export_directory=None, final=False):
        """
        Base export class. Overloading should really be avoided.

        Inputs:
            - naming_method: specifies how the output name is to be generated
                * 'raw_file'   = use the name and directory of the raw file, but with extension *.nc
                * 'last_load'  = use the name and directory of the last netcdf loaded
                * 'convention' = use the inbuilt naming convention for the filename. Directory will be last load if that is not blank, else it will be the raw_dir.
                * None         = use 'last_load' if it exists, else use 'raw_file'
            - export_directory: specify new directory if you don't want the program to choose
        """

        outdir, outname = self.get_export_location(naming_method=naming_method, export_directory=export_directory)

        self.attrs['last_export_directory'] = outdir
        self.attrs['last_export_file_name'] = outname

        ## 
        if final:
            self.ds.attrs['Disclaimer'] = self.disclaimer
            # outname = outname + 'finalised'

        self.ds.close() # Force close

        nc_file = '{folder}//{file}'.format(folder=outdir, file=outname)

        print('Exporting {}'.format(nc_file))
        
        self.ds.to_netcdf(path=nc_file)

        return self.ds # It may be useful in many cases to return somethng here. Just return self.sd for consistency with subclasses.  

    @property
    def best_export_dir(self):
        """
        Directory will be last load if that is not blank, else it will be the raw_dir.
        """

        outdir  = self.attrs['last_load_directory']
        if outdir == '':
            print('Best export dir is the raw_file_directory.')
            outdir  = self.attrs['raw_file_directory']
        else:
            print('Best export dir is the last_load_directory.')

        
        return outdir

    def get_export_location(self, naming_method=None, export_directory=None):
        """
        Base export class. Overloading should really be avoided.

        Inputs:
            - naming_method: specifies how the output name is to be generated
                * 'raw_file'   = use the name and directory of the raw file, but with extension *.nc
                * 'last_load'  = use the name and directory of the last netcdf loaded
                * 'convention' = use the inbuilt naming convention for the filename. Directory will be last load if that is not blank, else it will be the raw_dir.
                * None         = use 'last_load' if it exists, else use 'raw_file'
            - export_directory: specify new directory if you don't want the program to choose
        """

        if naming_method is None:
            if not len(self.attrs['last_load_directory']) == 0:
                naming_method = 'last_load'
            else:
                naming_method = 'raw_file'

        if naming_method.lower() == 'convention':
            print('Generating filename from naming convention.')
            outname = self.generate_outname()
            outname = zfile.drop_extension(outname) + '.nc'
            outdir  = self.best_export_dir

        elif naming_method.lower() == 'last_load':
            outdir  = self.attrs['last_load_directory']
            outname = self.attrs['last_load_file_name']
            print('Will save file using last load file name and directory.')

        elif naming_method.lower() == 'raw_file':
            outdir  = self.attrs['raw_file_directory']
            outname = self.attrs['raw_file_name']
            # Replace extension
            outname = zfile.drop_extension(outname) + '.nc'
            print('Will save file using raw file name and directory.')

        ## Use new directory if specified
        if not export_directory is None:
            outdir = export_directory
            print('Will to user input directory.')

        return outdir, outname

    @property
    def fullpath_last_export(self):
        """
        Full path that you last exported to.
        """

        return os.path.join(self.attrs['last_export_directory'], self.attrs['last_export_file_name'])

    def time_trim(self, first_good, last_good):

        if first_good is None:
            first_good = self.ds.time.values[0]

        if last_good is None:
            last_good = self.ds.time.values[-1]

        self.first_good = first_good
        self.last_good = last_good

        print('Good data from {} to {}'.format(first_good, last_good))

        # This sometimes fails if the DataSet time is timezone aware and the trim times
        # are timezone naive. Dolfyn seems to throw both:
            #  type(rr.ds.time.values[1]) = pandas._libs.tslibs.timestamps.Timestamp
            #  type(rr.ds.time.values[1]) = numpy.datetime64
        # and only the former is datetime aware. Ultimately I'll be going away from Dolfyn. 
        # there is no fix for this, user must be aware of this. 
        self.ds = self.ds.sel(time=slice(self.first_good, self.last_good))

        print('Trimmed Time')

    @property  
    def is_profile_data(self):
        
        if 'is_profile_data' in self.attrs:
            return self.attrs['is_profile_data']
        else:
            return False

    @property  
    def qc_conv(self):
        return self._qc_conv

    @qc_conv.setter
    def qc_conv(self, new_qc_conv):
        """
        Set function for _qc_conv
        """
        if isinstance(new_qc_conv, qc_conventions):
            self._qc_conv = new_qc_conv
        else:
            raise(Exception(("Must ba a QC Convention object")))

    def run_qc_compliance(self, flag_name):
        """
        Wrapper for the qc_conv.run_compliance method.

        Basically it checks that your QC flag [ds data_var with the name flag_name] complies with the QC conventions specified by the _qc_conv attribute.
        """

        self.qc_conv.run_compliance(self.ds[flag_name])

    def associate_qc_flag(self, var_name, flag_name):
        """
        Assign a QC flag to a DataArray, and add the flag to the dataset if necessary.
        Inputs:
                - var_name: The name of the data variable which will have a QC flag associated. Cannot QC dimensions [coords].
                - flag_name: The name of QC flag to be associated/added. The prefix 'qc_' will be added automatically
        """

        # Add the qc_ prefix
        flag_name = 'qc_' + flag_name
        
        # Check if var_name exists
        if not var_name in self.ds.data_vars:
            raise(Exception('Variable does not exist'))

        var_shape = self.ds[var_name].values.shape

        # Check if flag_name exists
        if not flag_name in self.ds.data_vars:
            # if not create the flag
            flag_array = self.ds[var_name].copy()
            flag_array.attrs={}
            flag_array.values[:] = -999
            self.ds[flag_name] = flag_array

            # flag = -999*np.ones_like(self.ds[var_name])
            # self.ds[flag_name] = flag
            self.qc_conv.run_compliance(self.ds[flag_name])
        elif not self.ds[flag_name].values.shape == var_shape:
            # If so check the size of the variable against the size of the flag
            raise(Exception("The selected QC flag does not have the right dimensions for this variable."))

        # Check this variable doesn't already have a qc code. Multiples not yet allowed.
        if 'qc_variable' in self.ds[var_name].attrs:
            raise(Exception("A QC variable has already been associated to the variable {}!".format(var_name)))

        for i in self.ds.data_vars:
            if 'is_qc_flag' in self.ds[i].attrs and 'associated data variables' in self.ds[i].attrs: # It might be a QC flag
                if  self.ds[i].attrs['is_qc_flag']: # It's a QC flag
                    if var_name in self.ds[i].attrs['associated_data_variables'].split(';'):
                        raise(Exception("The variable {} has already been associated to the QC flag {}!".format(var_name, i)))

        # Associate this flag to this variable
        self.ds[flag_name].attrs['associated_data_variables'] += ';'+var_name
        self.ds[var_name].attrs['qc_variable'] = flag_name

        # Trim leading ; if necessary
        if self.ds[flag_name].attrs['associated_data_variables'][0] == ';':
            self.ds[flag_name].attrs['associated_data_variables'] = self.ds[flag_name].attrs['associated_data_variables'][1::]
    
    @property  
    def valid_qc_flags(self):
        """
        All valid QC flags. Checks for:
            - 'qc_' prefix
            - presence of 'is_qc_flag' attribute
            - value of 'is_qc_flag' attribute

        The checks could certainly be more thorough.
        """

        out = []
        for i in self.ds.data_vars:
            if len(i) < 4:
                continue
            if not i[0:3] == 'qc_':
                continue
            if not 'is_qc_flag' in self.ds[i].attrs:
                continue
            if self.ds[i].attrs['is_qc_flag']:
                out.append(i)

        return(out)

    def set_pressure_sensor_height(self, ps_hasb):

        self.update_attribute('pressure_sensor_height_asb', ps_hasb)

    def get_pressure_sensor_height(self):

        pressure_sensor_height_asb = self._attrs['pressure_sensor_height_asb']

        if type(pressure_sensor_height_asb) == str:
            pressure_sensor_height_asb = np.nan
            
        return pressure_sensor_height_asb

    @property
    def pressure_sensor_height(self):

        return self.get_pressure_sensor_height()

    def advance_time_mins(self, advance_mins=None, dataset_name='ds', time_name='time', comment=""):
        """
        Function to shift clock and log the change. 

        Inputs:
            advance_mins is the number of minutes to advance the clock [i.e. move to a later time].
            dataset_name is the name of the dataset this change will be applied to. Default is "ds".
            time_name is the name of the time variable this change will be applied to. Default is "time".
             
        """

        if advance_mins == 0 or advance_mins is None or np.isnan(advance_mins):
            print('Not advancing time')
            return
            
        string = 'Advanced the time variable "{}" by {} minutes with user comment "{}"'.format(time_name, advance_mins, comment)
        self.add_comment(self.default_user, string, ds_name='ds', data_var=None)
        
        print(string)

        dataset = getattr(self, dataset_name)
        dataset = dataset.assign_coords({time_name: dataset[time_name] + np.timedelta64(advance_mins,'m')})
        setattr(self, dataset_name, dataset)

    def offset_variable(self, offset_in_var_units, variable_name, dataset_name='ds',  comment=""):
        """
        Modify raw data by adding a set offset. 

        The offset is sepcified in whatever units the variable is specified in.  

        """
    
        if offset_in_var_units == 0 or offset_in_var_units is None or np.isnan(offset_in_var_units):
            print('Not offsetting variable')
            return
        
        raise(Exception("Not implemented"))
        
    def overwrite_variable(self, new_val_in_var_units, variable_name, dataset_name='ds',  comment=""):
        """
        Modify variable by overwriting its value(s) with a constant or equal sized array.  
        
        The new value(s) is sepcified in whatever units the variable is specified in. 
        """
        
        string = 'Setting the variable "{}" to {} [same units] with user comment "{}"'.format(variable_name, 
                                                                                              new_val_in_var_units, comment)
        self.add_comment(self.default_user, string, ds_name='ds', data_var=None)
        
        print(string)

        dataset = getattr(self, dataset_name)
        dataset[variable_name].values[:] = new_val_in_var_units
        

    def update_qc_flag_dict(self, flag_name, index_dict, flag_value, comment=None, delete_raw=False, verbose=False):
        """
        This is a base function to update a QAQC flag. Inputs:
            - flag_name: The name of the netcdf variable corrsponding to the QC flag being edited. Can be:
                            - a single flag name
                            - a list of flag names
                            - '*'
                        if '*' is used all QC flags will be updated 

            - index_name: The name of the netcdf variable which is being used to as an index to 
                            identify which points to in the QC flag variable are to be edited.
            - start: The first point in the index to change the value of
            - end: The last point in the index to change the value of
            - flag_value: new value which the QC flag is to become wherever the index is between start and end 
        """
        
        if delete_raw:
            error

        if isinstance(flag_name, str): 
            if flag_name == '*': 
                flag_names = self.valid_qc_flags
            else:
                flag_names = [flag_name]
        else: # Assume list
            flag_names = flag_name

        for flag_name in flag_names:

            flag_dims = [dim for dim in self.ds[flag_name].dims]
            index_check = [index_name in self.ds[flag_name].dims for index_name in index_dict.keys()]
            if not np.all(index_check):
                continue

            # logind = np.ones_like(self.ds[flag_name].values.astype(bool))
            # for index_name in index_dict.keys():

            #     index = expand_dims(self.ds[index_name], flag_dims)

            #     start, end = index_dict[index_name]

            #     logind1 = index >= start
            #     logind = np.logical_and(logind, logind1)

            #     logind2 = index <= end
            #     logind = np.logical_and(logind, logind2)

            # # Note that the logical index is reversed here.
            # self.ds[flag_name] = self.ds[flag_name].where(~logind, flag_value)
                
            my_count = self.ds[flag_name].copy()
            my_count.values[:] = 0
            conditions = 0
            for index_name, index_vals in index_dict.items():
                
                start, end = index_vals
                
                mask_0D = (
                    (self.ds.coords[index_name] >= start) &
                    (self.ds.coords[index_name] <= end) 
                )

                print('SUM OF MASK = {}'.format(np.sum(mask_0D.values)))
                
                if False:
                    my_count[mask_0D] += 1 # Doesn't work for 2D
                else:
                    dims = self.ds[flag_name].dims    # Get a list of dims
                    hard = [np.s_[:] for dim in dims] # Make a list with full slice as default
                    mydim = [i for i in np.arange(len(dims)) if dims[i] == index_name][0] # Find our index axis
                    hard[mydim] = mask_0D.values # Set this axis to the mask
                    hard = tuple(hard) # Convert to tuple because it just works
                    
                    my_count[hard] += 1 # Incriment
                    
                # my_count.sel({index_name: mask_0D}) # Doesn't work at all. Ask Professor xarray about this one.  

                conditions += 1
                                
            logind = my_count.values == conditions
            logind = ~logind
            
            print('SUM OF logind = {}'.format(np.sum(logind)))
            print('SUM OF ~logind = {}'.format(np.sum(~logind)))

            self.ds[flag_name] = self.ds[flag_name].where(logind, flag_value)
            
            ind = np.where(logind)[0]
            logind = self.ds[flag_name] > 0

            if not comment is None: # Log comment on the QAQC.
                string = 'Flagged {} values with code "{}" and user comment "{}"'.format(len(self.ds[index_name].values)-len(ind), flag_value, comment)

                self.add_comment(self.default_user, string, data_var=flag_name)

            if self.verbose:
                
                print(ind)
                print(len(ind))
                
                print('There are {} flagged data points.'.format(np.sum(logind.values)))
                print('')        

    def update_qc_flag(self, flag_name, index_name, start, end, flag_value, comment=None, delete_raw=False, verbose=False, allow_back_compat=True):
        """
        This is a base function to update a QAQC flag. Inputs:
            - flag_name: The name of the netcdf variable corrsponding to the QC flag being edited. Can be:
                            - a single flag name
                            - a list of flag names
                            - '*'
                        if '*' is used all QC flags will be updated 

            - index_name: The name of the netcdf variable which is being used to as an index to 
                            identify which points to in the QC flag variable are to be edited.
            - start: The first point in the index to change the value of
            - end: The last point in the index to change the value of
            - flag_value: new value which the QC flag is to become wherever the index is between start and end 
        """
        
        if not allow_back_compat:
            raise(Exception('The function update_qc_flag is no longer available. Use the function update_qc_flag_dict or try again with the allow_back_compat option. '))
            
        index_dict = {}
        index_dict[index_name] = [start, end]

        self.update_qc_flag_dict(flag_name, index_dict, flag_value, comment=comment, delete_raw=delete_raw, verbose=verbose)

    def update_qc_flag_logical(self, flag_name, index_name, logical_index, flag_value, comment=None, delete_raw=False, verbose=False):
        """
        This is a base function to update a QAQC flag. Inputs:
            - flag_name: The name of the netcdf variable corrsponding to the QC flag being edited. Can be:
                            - a single flag name
                            - a list of flag names
                            - '*'
                        if '*' is used all QC flags will be updated 

            - index_name: The name of the netcdf variable which is being used to as an index to 
                            identify which points to in the QC flag variable are to be edited.
            - logical_index: Logical index for points to change the value of
            - flag_value: new value which the QC flag is to become wherever the index is between start and end 
        """

        if delete_raw:
            error

        if isinstance(flag_name, str): 
            if flag_name == '*': 
                flag_names = self.valid_qc_flags
            else:
                flag_names = [flag_name]
        else: # Assume list
            flag_names = flag_name

        for flag_name in flag_names:

            if not index_name in self.ds[flag_name].coords:
                continue

            # Note that the logical index is reversed here.
            self.ds[flag_name] = self.ds[flag_name].where(~logical_index, flag_value)

            ind = np.where(logical_index)[0]
            logind = self.ds[flag_name] > 0

            if self.verbose:
                print()

            if not comment is None: # Log comment on the QAQC.
                # string = 'Flagged {} values with code "{}" and user comment "{}"'.format(len(self.ds[index_name].values)-len(ind), flag_value, comment)
                string = 'Flagged {} values with code "{}" and user comment "{}"'.format(len(ind), flag_value, comment)
                    
                self.add_comment(self.default_user, string, data_var=flag_name)

            if self.verbose:
                print(ind)
                print(len(ind))

                print('There are {} flagged data points.'.format(np.sum(logind.values)))
                print('')   

    def flip_qc_value(self, flag_name, value_out=-999, value_in=0, verbose=False):
        """
        Flip one QC value for another. By default it will flip -999 to 0

        flag_name can be:
                                - a single flag name
                                - a list of flag names
                                - '*'
                            if '*' is used all QC flags will be updated 
        
        NOTE: doesn't actually check that the input variable is QC convention compliant. 
        
        """
        
        if isinstance(flag_name, str): 
                if flag_name == '*': 
                    flag_names = self.valid_qc_flags
                else:
                    flag_names = [flag_name]
        else: # Assume list
            flag_names = flag_name

        for flag_name in flag_names:

            logind = self.ds[flag_name] == value_out
            logind = ~logind

            self.ds[flag_name] = self.ds[flag_name].where(logind, value_in)

            if self.verbose:
                self.get_fig_text()

                logind = self.ds[flag_name] > 0
                print('There are {} positive flag points.'.format(np.sum(logind.values)))
                logind = self.ds[flag_name] < 0
                print('There are {} negative flag points.'.format(np.sum(logind.values)))
                print('')   

    def get_qaqc_var(self, var_name):
        """
        Retrun a QAQC'd copy of the data array by var_name. This just routes to the main function.
        """
        
        da = get_qaqc_var(self.ds, var_name)
        
        return da 

    def has_dates(self):

        if type(self.so) == type(None):
            return False
        elif type(self.eo) == type(None):
            return False
        else:
            return True

    def get_fig_text(self, fileappend=''):
        """
        Get a text string to describing the dataset.
        """

        delim = ' | '
        s = ''
        k = self.attrs.keys()
        
        if 'project' in k:
            s = s + delim + self.attrs['project']
            
        if 'trip_recovered' in k:
            s = s + delim + self.attrs['trip_recovered']
            
        if 'site' in k:
            s = s + delim + self.attrs['site']
            
        if not fileappend == '':
            s = s + delim + fileappend
            
        s = s[len(delim)::]

        print(s)

        return s

    def add_history(self, author, string, ds_name=None, data_var=None):
        """
        Add CF Compliant string to the history attribute of a Dataset or DataArray.

        Use the data_var kwarg to specify a DataArray. Omit or set to None to work with the whole Dataset.
        """

        attr = 'history'
        self.add_string(attr, author, string, ds_name=ds_name, data_var=data_var)
            
    def add_comment(self, author, string, ds_name='ds', data_var=None):
        """
        Add CF Compliant string to the comments attribute of a Dataset or DataArray.

        Use the data_var kwarg to specify a DataArray. Omit or set to None to work with the whole Dataset.
        """

        attr = 'comment'
        self.add_string(attr, author, string, ds_name=ds_name, data_var=data_var)
            
    def add_string(self, attr, author, string, ds_name='ds', data_var=None):
        """
        Add CF Compliant string to an attribute of a Dataset or DataArray.

        Use the data_var kwarg to specify a DataArray. Omit or set to None to work with the whole Dataset.
        """
        
        if data_var is None:
            obj = getattr(self, ds_name)
        else:
            obj = getattr(self, ds_name)[data_var]

        new_string = datetime.datetime.now().isoformat() + ': ' + '[{}]'.format(author) + ' ' + string
        
        if not attr in obj.attrs.keys():
            obj.attrs[attr] = new_string
        elif obj.attrs[attr] in ['', '?', 'blank']:
            obj.attrs[attr] = new_string
        else:
            obj.attrs[attr] = obj.attrs[attr] + ';' + new_string
            
        return new_string

    def update_attribute(self, attribute_name, attribute_value, ds_name='ds', data_var=None, strict=True):
        """
        This function updates the hidden attributes property of the class. The attribute must exist in the default attributes dictionary.
        """

        if data_var is None:
            obj = getattr(self, ds_name)
        else:
            obj = getattr(self, ds_name)[data_var]

        if (attribute_name in self._default_attrs) or (not strict):
            if self.verbose:
                print('Setting attribute "{}" to "{}"'.format(attribute_name, attribute_value))
            obj.attrs[attribute_name] = attribute_value
        else:
            raise(Exception('{} is not a valid attribute.'.format(attribute_name)))

        # Repeated calls are made here, shouldn't be too slow. 
        self.pIMOS_assign_coords()

    def update_attributes_with_dict(self, attribute_dict):
        """
        This function updates the hidden attributes property of the class. The attribute must exist in the default attributes dictionary.
        """

        for attribute_name in attribute_dict.keys():
            
            self.update_attribute(attribute_name, attribute_dict[attribute_name])

    
    def pIMOS_assign_coords(self):
        """
        Assign some additional coords to the file. These coords use deployment knowledge, so are not available at first creation of the ds.  

        Note: at the moment this doesn't alter the coords on any data variables.

        Function can be called repeatedly, it simply overwrites any coords.

        """
            
        try:
            lat_nom = self.ds.attrs['nominal_latitude']
            self.ds = self.ds.assign_coords({'lat_nom': lat_nom,})
        except:
            # print('Failed assigning nominal_latitude')
            pass
            
        try:
            lat_nom = self.ds.attrs['nominal_longitude']
            self.ds = self.ds.assign_coords({'lon_nom': lat_nom,})
        except:
            # print('Failed assigning nominal_longitude')
            pass

        try:
            z_nom = self.ds.attrs['nominal_site_depth'] + self.ds.attrs['nominal_instrument_height_asb'] 
            self.ds = self.ds.assign_coords({'z_nom': z_nom,})
        except:
            # print('Failed assigning z_nom')
            pass
            
    
    @property
    def rawattrs(self):
        """
        Return the raw file attribute string as a dict.
        
        This is sloppy, and will likely fail under many data types. 
        """
        
        array = np.array
        
        self.attrs['raw_file_attributes'] = "{'_type': 'ADCP', 'name': 'unrecognized firmware version', 'sourceprog': 'instrument', 'prog_ver': 51.4, 'config': '01000001-11001010', 'beam_angle': 20, 'numbeams': 4, 'beam_freq_khz': 300, 'beam_pattern': 'convex', 'orientation': 'up', 'simflag': 'real', 'n_beam': 4, 'n_cells': 46, 'pings_per_ensemble': 60, 'cell_size_m': 2.0, 'blank_m': 1.76, 'prof_mode': 1, 'corr_threshold': 0, 'prof_codereps': 5, 'min_pgood': 0, 'evel_threshold': 2000, 'sec_between_ping_groups': 1.0, 'coord': '00000000', 'coord_sys': 'beam', 'use_pitchroll': 'no', 'use_3beam': 'no', 'bin_mapping': 'no', 'xducer_misalign_deg': 0.0, 'magnetic_var_deg': 0.0, 'sensors_src': '01111101', 'sensors_avail': '00111101', 'bin1_dist_m': 4.2, 'xmit_pulse': 2.38, 'water_ref_cells': array([1, 5]), 'fls_target_threshold': 255, 'xmit_lag_m': 0.49}"
        
        return eval(self.attrs['raw_file_attributes'])
        
    def get_raw_file_attributes(self):
        """
        Returns the attributes from rawfile as a dict.
        """

        return eval(self.attrs['raw_file_attributes'])

    def store_raw_file_attributes(self, ds):
        """
        Function that will pull any valid attributes from the dataset into the wrapper, and will
        push any remaining into the 'raw_file_attributes' attribute.
        """
        
        print('STORING RAW FILE ATTRIBUTES')
        input_attributes = ds.attrs.copy()
        keys = input_attributes.keys()

        invalid_attributes = {} # These can't be added to the dataset 
        ds.attrs = self._default_attrs  
        
        for attribute in keys:
            
            if attribute in xrwrap._default_attrs:
                self.update_attribute(attribute, input_attributes[attribute])
            else:
                invalid_attribute = input_attributes[attribute]
                invalid_attributes[attribute] = invalid_attribute
                print(attribute.upper())
        pass

        if not len(invalid_attributes) == 0:

            if not ds.attrs['raw_file_attributes'] == '':
                raise(Exception('This operation can be performed only once'))

            # self.add_string('raw_file_attributes', self.default_user, str(invalid_attributes), ds_name='ds', data_var=None)
            self.attrs['raw_file_attributes'] = str(invalid_attributes)

    def enforce_these_attrs(self, enforced_attributes):
        """
        Enforces certain global attributes to have certain values. Error will be thrown if the attribute are not
        blank and different to what is specified:

        Inputs:
            - enforced_attributes (dict): dictionary of with attribute names as keys, and the enforced value as value.  
        """

        # CF Compliance
        if self.verbose:
            print(self.attrs)

        for attr in enforced_attributes.keys():
            if self.attrs[attr] == '':
                self.update_attribute(attr, enforced_attributes[attr])
            else:
                assert(self.attrs[attr]==enforced_attributes[attr])  

    def check_attrs(self):
        """
        Check that all attributes in the dataset are in the default attributes, and vice versa.
        
        TO DO:
            - There should probably be a list of attributes that need not be blank or a warning given. 
        """
        
        for key in self.ds.attrs.keys():
            assert(key in self._default_attrs)
        for key in self._default_attrs:
            assert(key in self.ds.attrs.keys())
    
    def parse_attributes(self, ds_to_check=None):
        """
        Generic function to check consistency of attributes between the object itself and the properties of the class.
        """

        raise(Exception('We no longer store attributes this way - this function should not be used'))

        if ds_to_check is None:
            ds_to_check = self.ds

        print("Parsing attributes.")
        for i in ds_to_check.attrs.keys():
            if i in self.attrs.keys():
                print("{} is both a property of the object and an attribute of the dataset".format(i))
                if ds_to_check.attrs[i] == self.attrs[i]:
                    print("     ... and they are equal")
                else:
                    print("     ... and they NOT are equal!!!")

        ds_to_check.attrs = self.attrs

    @property
    def disclaimer(self):
        """

        """

        disc = """These data were prepared for a specific purpose at a specific site, which may or may not be disclosed in the data or metadata. 
        Use of these data for any purpose, without the express written consent of the author is not permitted. 
        The data are provided WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
        The author is not liable in any way for consequences arising from any use of these data or any derivatives in any application including but not limited to design or decision-making processes.
        """

        return disc
  
def get_qaqc_var(ds, var_name):
    """
    Retrun a QAQC'd copy of the data array by var_name
    """
    
    da = ds[var_name].copy() # return a QCd copy

    if not 'qc_variable' in ds[var_name].attrs:
        print('Variable {} has no QAQC flag. Returning raw.'.format(var_name))
    else:    
        flag_name = ds[var_name].attrs['qc_variable'] # Currently assumes only one qc_variable
        ind = ds[flag_name].values > 0 # Currently does not remove suspect data. Should be an option.
        print('Blanking {} values.'.format(np.sum(ind)))

        values=da.values
        values[ind] = np.nan
        da.values=values
    
    return da             

### Generic XR stuff

def calendar_month(year, month):
    """
    For a given year and month return the date of the begining of the month and the date of the beginning of the next month
    """
    start = datetime.datetime(year, month, 1)
    if month == 12:
        end = datetime.datetime(year+1, 1, 1)
    else:
        end = datetime.datetime(year, month+1, 1)
    print(start)
    print(end)
    return start, end
    
def select_calendar_month(X, year_month, timename='time'):
    """
    Return a subset of the Dataset or DataArray limited to the specified year and month. Uses rather crude indexing methods. I know 
    there are much more sopistocated inbuilt functions but I'm not familiar with them.
    """

    
        
    year, month = year_month
    
    start, end = calendar_month(year, month)
    
    # Gotta do better than this
    if timename.lower() == 'time':
        X_cm = X.sel(time = slice(start, end))
    elif timename.lower() == 'time_wave':
        X_cm = X.sel(time_wave = slice(start, end))
    # Gotta do better than this
    
    return X_cm, [start, end]

def split_by_timegap(X, timename='time', hours=1):
    """
    Split a DataArray along the time dimension based on a set time gap.
    
        Returns a list of DataArrays with the gaps removed.
        
    """
        
    time = X[timename].values
    dt = np.diff(time)

    # print(len(time))
        
    i = [int(ii) for ii in np.where(dt>np.timedelta64(hours, 'h'))[0]] + [len(time)-1]
    i = [-1] + list(set(i))
    i.sort()
    
    # print('Split index')
    # print(i)
    
    Xs = [X.isel({timename: np.arange(i[j]+1, i[j+1])}) for j in np.arange(0, len(i)-1)]
    
    return Xs

# def expand_dims(array, full_dims):
#     """
#     Expands out an array to help with broadcasting.

#     Inputs: 
#         - array: Array you wish to expand
#         - full_dims: List containing names of dimensions of the array you wish to expand out to

#     Returns:
#         - vals: values of the array expanded into the dimensions within full_dims

#     """
    
#     expand_list = []
#     for flag_dim in full_dims:
#         if flag_dim in array.coords:
#             expand_list.append(np.arange(0, len(array[flag_dim])))
#         else:
#             expand_list.append(None)

#     vals = array.values[expand_list]

    return vals

def get_dim_index(array, dim_name):
    """
    Find which axis of a DataArray corresponds to a certain dimension name.
    """
        
    coord_names =[x for x in array.coords ]
    dim_index =[i for i in np.arange(0, len(coord_names)) if coord_names[i]==dim_name]    
    
    if len(dim_index) == 1:
        dim_index=dim_index[0]
    else:
        dim_index = None
        
    return dim_index

def naming_conv(attrs, convention=None):
    """
    Create a 'conventional name' from a list of attributes defining the naming convention.      
    """
    
    if convention is None:
        convention = ['project',
                    'trip',
                    'site_station',
                    'instrument_model',
                    'instrument_serial_number']
    
    name = ''
    for attribute in convention:
        val = attrs[attribute]
        if type(attrs[attribute]) == str:
            # Check if it has been ; delimited, in which case return a special output
            split_attr = val.split(';')
            split_attr = list(set(split_attr)) # Just the unique ones
            n = len(split_attr)
            if n == 1:
                val = split_attr[0]
            else:
                val = 'multi({})_'.format(n)+attribute

        elif np.isnan(val):
            val = ''
        else:
            val = str(int(attrs[attribute]*100))
            
        if val == '':
            val = 'no_'+attribute

        if '.' in val:
            val = zfile.drop_extension(val)

        name += '[{}]_'.format(val)
    
    name = name[0:-1]
    
    return name


# %%
