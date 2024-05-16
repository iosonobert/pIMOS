import pIMOS
import os, glob
import pandas as pd, numpy as np, xarray as xr



def get_pimos_catalogue(root_folder, is_moored=True):
    """
    Loads the PIMOS catalogue from a specified root folder. The root folder should contain a subfolder
    called 'pimos_vX.Y', where X.Y is the version of the pIMOS library you are using. 
        
    The function iterates over all subfolders and files in the given root folder, 
    checking for compliance with a specific folder structure and file naming convention. 
    It extracts metadata from compliant file names and stores it in a pandas DataFrame.

    **HAVING WRITTEN THIS FUNCTION, I SEE WE COULD HAVE DONE A BETTER JOB WITH THE NAMING CONVENTIONS AND FOLDER SETUP**

    Parameters
    ----------
    root_folder : str
        The path to the root folder containing the PIMOS catalogue.
    is_moored : bool, optional
        If True, the function assumes that the data is from moored instruments and 
        does not output the raw filename. If False, the function assumes that the data 
        is from profiling instruments and does not output the nominal height above bed. 
        Default is True. **IMPORTANT BECAUSE NAMING CONVENTION DIFFERS**

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the PIMOS catalogue. The columns of the DataFrame are:
        'File version', 'Group', 'project', 'trip', 'site_station', 
        'nominal_instrument_height_asb' (moored only), 'instrument_model', 
        'instrument_serial_number', 'raw_file_name' (profiling only), 'outfile_append', 'nc_path'.

    Raises
    ------
    AssertionError
        If the PIMOS folder does not exist in the root folder.
    Warning
        If a folder or file does not comply with the expected structure or naming convention, 
        a warning is printed to the console.
    """

    version_subfolder = f'pimos_v{pIMOS.__version__}'

    print(f'Looking for {version_subfolder} within {root_folder}')

    pimos_folder = os.path.join(root_folder, version_subfolder)

    assert os.path.exists(pimos_folder), f'Folder {pimos_folder} does not exist'

    print(f'.... {version_subfolder} found')

    folders = glob.glob(pimos_folder+'/*')
    fvfolders = os.listdir(pimos_folder)

    default_row = np.array(['N/A']*8).astype('<U30')
    catalogue = []

    # Ideally would load this from a module somewhere. For a later version perhaps.
    columns = ['project',
                'trip',
                'site_station',
                'nominal_instrument_height_asb',   # Moored only
                'instrument_model',
                'instrument_serial_number',
                'raw_file_name',                   # Profiling only
                'outfile_append']

    all_columns = ['File version', 'Group'] + columns + ['nc_path']

    for fvfolder in fvfolders:
        fvpath = os.path.join(pimos_folder, fvfolder)
        print(fvfolder)

        # Couple of base checks
        if not len(fvfolder)>2:
            print(f'WARNING: folder structure non compliant {fvpath}')
        if not fvfolder.lower()[0:2]=='fv':
            print(f'WARNING: folder structure non compliant {fvpath}')

        # groupsfolders = glob.glob(fvpath+'/*')
        groupsfolders = os.listdir(fvpath)

        for groupfolder in groupsfolders:
            # print(groupfolder)
            groupfolderpath = os.path.join(fvpath, groupfolder)
            
            nc_paths = glob.glob(groupfolderpath+'/*.nc')
            # nc_paths = os.listdir(groupfolderpath)
            # print(nc_paths)
            for nc_path in nc_paths:
                nc_folder, nc_file = os.path.split(nc_path)
                nc_file = os.path.splitext(nc_file)[0]

                compliant = True
                if not len(fvfolder)>2:
                    compliant = False
                if not nc_file[0] == '[':
                    compliant = False
                if not nc_file[-1] == ']':
                    compliant = False

                if not compliant:
                    print(f'WARNING: folder structure non compliant {nc_path}')
                    continue
                
                nc_file = nc_file[1:-1]
                nc_file_parts = nc_file.split(']_[')
                # print(nc_file_parts)

                my_row = default_row.copy()
                # Now the annoying part is we have no way to confirm the naming convention used
                if is_moored:
                    # Won't output column 6 - the raw filename
                    my_row[[0, 1, 2, 3, 4, 5, 7]] = nc_file_parts
                else:
                    # Won't output column 3 - the nominal height above bed
                    my_row[[0, 1, 2, 4, 5, 6, 7]] = nc_file_parts
                

                full_row = [fvfolder, groupfolder] +  list(my_row) + [nc_path]

                catalogue.append(full_row)

    return pd.DataFrame(catalogue, columns=all_columns) 


def filter_catalogue(cat, filters, reverse=False):
    """
    Filters a catalogue based on a dictionary of filters.

    Parameters
    ----------
    cat : pandas.DataFrame
        The DataFrame to filter.
    filters : dict
        A dictionary where the keys are column names and the values are the values to filter by.
        The values can be either single values or lists of values.

    Returns
    -------
    cat_f = pandas.DataFrame
        The filtered DataFrame.
    """

    cat_f = cat.copy()
    for key, value in filters.items():
        if isinstance(value, list):
            if reverse:
                cat_f = cat_f[~cat_f[key].isin(value)]
            else:
                cat_f = cat_f[cat_f[key].isin(value)]
        else:
            if reverse:
                cat_f = cat_f[cat_f[key] != value]
            else:
                cat_f = cat_f[cat_f[key] == value]

    return cat_f


def load_detailed_data(cat_f):
    """
    Adds detailed metadata to a preloaded a PIMOS catalogue DataFrame. The data added are within the 
    file itself, and so the file must be potentiall downloaded and opened. This is potentially a slow 
    operation, so it is recommended to do some filtering before calling this function.

    Metadata added are:
    - nominal latitude
    - nominal longitude
    - start time [total record, not first good data]
    - end time   [total record, not first good data]


    Parameters
    ----------
    cat_f : pandas.DataFrame
        The PIMOS catalogue DataFrame.

    Returns
    -------
    pandas.DataFrame
        The updated PIMOS catalogue DataFrame with additional columns for nominal latitude, 
        nominal longitude, start time, and end time.
    """

    nominal_latitude = []
    nominal_longitude = []
    start_time = []
    end_time = []
    for i, row in cat_f.iterrows():

        # if False: # this seems unnecessary
        #     rr = pIMOS.load_pimos_nc(row['nc_path'])
        #     ds = rr.ds
        #     attrs = rr.attrs
        # else:
        ds = xr.open_dataset(row['nc_path'])
        attrs = ds.attrs

        nominal_latitude.append(attrs['nominal_latitude'])
        nominal_longitude.append(attrs['nominal_longitude'])

        start_time.append(ds.time.values[0])
        end_time.append(ds.time.values[-1])


    cat_f['nominal_latitude'] = nominal_latitude
    cat_f['nominal_longitude'] = nominal_longitude
    cat_f['start_time']       = start_time
    cat_f['end_time']         = end_time

    return cat_f


def get_var_bounds(catalogue, var):
    """
    Get the bounds of a variable in a catalogue.

    Parameters
    ----------
    catalogue : pandas.DataFrame
        The PIMOS catalogue DataFrame.
    var : string
        The variable to get the bounds of.

    Returns
    -------
    tuple
        A tuple containing the minimum and maximum values of the variable.
    """
    bounds_col = []
    min_val = np.nan
    max_val = np.nan
    for i, row in catalogue.iterrows():
        ds = xr.open_dataset(row['nc_path'])
        if (var in ds.data_vars) | (var in ds.coords):
            min_val = np.nanmin(ds[var].values)
            max_val = np.nanmax(ds[var].values)
        else:
            if ('nominal_'+var) in ds.attrs:
                min_val = np.nanmin(ds.attrs['nominal_'+var])
                max_val = np.nanmax(ds.attrs['nominal_'+var])
        bounds_col.append([min_val, max_val])
    catalogue[(var + '_bounds')] = bounds_col
    return catalogue


def ds_check(catalogue, var, ds_key, strict=True, filter_name=None):
    '''
    Read coordinates from each row in the catalogue and return boolean mask for each row if coordinate exists
    
    Parameters
    ----------
    catalogue : pandas.DataFrame
        The PIMOS catalogue DataFrame.
    coord : string
        The coordinate to filter by.
    strict : bool, optional
        If True, the function will return True if a partial match is found. E.g. 'vel' is in 'velocity'.
    filter_name : string, optional
        The name of the filter to add to the catalogue DataFrame.

    Returns
    -------
    pandas.DataFrame
        The updated PIMOS catalogue DataFrame with an additional column(s) for the filter.
    '''

    if filter_name is None:
        filter_name = f'{var}_exists'

    check_list = []

    # Loop through each row in the catalogue
    for i, row in catalogue.iterrows():

        # Load the netCDF file
        ds = xr.open_dataset(row['nc_path'])

        if ds_key == 'coords':
            variables = list(ds.coords)
        elif ds_key == 'variables':
            variables = list(ds.data_vars)
        elif ds_key == 'attrs':
            variables = list(ds.attrs)
        else:
            raise ValueError(f'Invalid ds_key: {ds_key}')

        # Check if coordinate is in the dataset coords
        check_list.append(dataset_key_check(variables, var, strict=strict))

    catalogue[filter_name] = check_list
    return catalogue



def get_data_cube(catalogue, cube_bounds, columns=['latitude_bounds', 'longitude_bounds', 'time_bounds']):
    '''
    Get a data cube from the catalogue based on the bounds of the coordinates
    '''
    check_list = True
    # Find if data is within the cube bounds
    assert len(cube_bounds) == len(columns), 'Number of cube bounds must match number of columns'
    for cub, col in zip(cube_bounds, columns):
        # Check if the coordinate is within the bounds
        check = (cub[0] <= catalogue[col].apply(lambda x: x[1])) & (cub[1] >= catalogue[col].apply(lambda x: x[0]))
        check_list = check_list & check
    return catalogue[check_list]


def dataset_key_check(key_list, key, strict=True):

    # Check if key is in the ist of keys
    if strict:
        if key in key_list:
            return True
        else:
            return False
    else:
        if np.any([key in v for v in key_list]):
            return [v for v in key_list if key in v]
        else:
            return ''