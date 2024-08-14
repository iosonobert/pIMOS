import datetime
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
import afloat.plot.plotting as zplot
import warnings
from pIMOS import __version__

"""
Various tools that should be pIMOS tools, but that I have not officially found a spot for yet. 
"""

def parse_dbconfig(dbconfig_file):
    """
    This just reads all of the files paths out of the config file
    """

    with open(dbconfig_file) as f:
        lines = f.readlines()

    config = {}    
    dbconfig_file_split = os.path.split(dbconfig_file)
    config['db_root'] = dbconfig_file_split[0]

    for line in lines:
        if len(line)<3:
            continue
        elif line[0] == '#':
            continue
        elif not '=' in line:
            continue
        else:
            print(line)

        words = line.split('=')

        if not len(words)==2:
            continue

        config[words[0].strip()] = words[1].strip()
    
    return config

def read_db(dbconfig_file):
    
    db_config = parse_dbconfig(dbconfig_file)

    db_data = {}
    for log in ['autonomous_metadata', 'profiler_metadata']:
        if log in db_config:
            db_data[log] = parse_db_csv(db_config, table_name=log)
        else:
            raise(Exception('Config file incomplete'))
            
    if 'possible_mooring_dates' in db_config:
        db_data['possible_mooring_dates'] = parse_possible_mooring_dates(db_config, recovered=None)
    else:
        raise(Exception('Config file incomplete'))
        
    
        
        
    return db_data

def parse_db_csv(db_config, table_name):

    file = os.path.join(db_config['db_root'], db_config[table_name])
    df = pd.read_csv(file)
    
    return df

def parse_deployment_metadata(db_config):

    # file = os.path.join(db_config['db_root'], db_config['deployment_metadata'])
    # df = pd.read_csv(file)
    df = parse_db_csv(db_config, table_name='deployment_metadata')

    return df

def parse_possible_mooring_dates(db_config, recovered=None):

    possible_mooring_dates_file = os.path.join(db_config['db_root'], db_config['possible_mooring_dates'])
    df = pd.read_csv(possible_mooring_dates_file, parse_dates=['StartDate', 'EndDate'], dayfirst=True)
    
    if not recovered is None:
        df = df.loc[df['Recovered'] == recovered]

    return df

def strcmpi(lst, string):
    """
    Copy of the matlab function.
    Returns list.
    """
    rtn = [i.lower() == string.lower() for i in lst]
    return rtn
    pass

def plot_echo(rr, db_data, mooring, attributes, variable='echo', width=65,\
             cmap='magma', experiment=None, recovered=None):
    
    #%matplotlib inline
    fig = plt.figure(figsize=(20,3))

    zl = zplot.axis_layer(left=2, right=2, bottom=2, top=2, heights = [4], widths=[width])
    zl.verbose = False

    zl.lay(0, 0)

    rr.get_qaqc_var(variable).plot(label='QAQC')

    plt.xlim(rr.ds.time.values[[0, -1]])
    rr.get_qaqc_var(variable)[0, :, :].plot(label='QAQC', cmap=cmap)

    plt.grid()

    title= ' | '.join([str(attributes[i]) for i in attributes])
    plt.title(title)

    add_mooring_dates(db_data, mooring, plt.gca(), experiment=experiment, recovered=recovered)
    plt.show()

    return fig


def plot_temp(rr, db_data, mooring, attributes, variable='Temperature',\
              plotraw=True, width=65, experiment=None, recovered=None, transpose=False):
    
    if experiment is None:
        try:
            experiment = attributes['project']
        except:
            Exception("Please supply an experiment.")

    if recovered is None:
        try:
            recovered = attributes['trip']
        except:
            Exception("Please supply a recovery trip name.")            
    
    #%matplotlib inline
    fig = plt.figure(figsize=(20,3))

    zl = zplot.axis_layer(left=2, right=2, bottom=2, top=2, heights = [4], widths=[width])
    zl.verbose = False

    zl.lay(0, 0)

    # Plotting
    data = rr.ds[variable].T if transpose else rr.ds[variable]
    qaqc_data = rr.get_qaqc_var(variable).T if transpose else rr.get_qaqc_var(variable)
    if plotraw:
        data.plot(label='Raw')
    qaqc_data.plot(label='QAQC')

    plt.xlim(rr.ds.time.values[[0, -1]])
    if plotraw:
        plt.legend()
    
    plt.grid()

    title= ' | '.join([str(attributes[i]) for i in attributes])
    plt.title(title)
    
    if not rr.attrs['is_profile_data']:
        if db_data is not None:
            add_mooring_dates(db_data, mooring, plt.gca(), experiment=experiment, recovered=recovered)
        else:
            print('No mooring dates available')
    plt.show()
    
    return fig

def add_mooring_dates(db_data, mooring, ax, experiment=None, recovered=None):
    
    yl = ax.get_ylim()
    
    df = db_data['possible_mooring_dates']

    # Hard Code This For Now
    df = df.loc[strcmpi(df['Recovered'].values, recovered)]
    df = df.loc[strcmpi(df['Experiment'].values, experiment)]
    
    if df.shape[0] == 0:
        raise(Exception("Could not find dates for this mooring. "))

    df = df.loc[strcmpi(df['Mooring'].values, mooring)]

    ax.plot([df['StartDate'].values[0]]*2, yl, 'r--')
    ax.plot([df['EndDate'].values[0]]*2, yl, 'r--')
    
    pass

def pIMOS_get_archive_folder(archive_dir, file_version):
    """
    Get the specific archive folder for the pIMOS version, and file process_level. 
    """

    folder = os.path.join(archive_dir, 'pimos_v'+__version__)
    if not os.path.exists(folder):
        os.mkdir(folder)

    if type(file_version) is str:
        file_version = 0

    fv = 'FV{:02.0f}'.format(file_version)
    folder = os.path.join(folder, fv)
    if not os.path.exists(folder):
        os.mkdir(folder)

    return folder

def pIMOS_get_archive_folder_model(archive_dir, file_version, model):

    folder = pIMOS_get_archive_folder(archive_dir, file_version)
    
    folder = os.path.join(folder, model)
    if not os.path.exists(folder):
        os.mkdir(folder)

    return folder

def pIMOS_get_my_export_folder(rr, archive_dir):

    model        = rr.ds.attrs['pimos_nickname']
    file_version = rr.ds.attrs['process_level']
    folder       = pIMOS_get_archive_folder_model(archive_dir, file_version, model)

    return folder

# def pIMOS_export_old(rr, archive_dir, model, file_append=''):

#     file_version = rr.ds.attrs['process_level']
#     file_version = rr.ds.attrs['process_level']
#     folder = pIMOS_get_archive_folder_model(archive_dir, file_version, model)

#     print(folder)

#     # rr.folder = folder
#     # rr.file_ = serial
    
#     rr.export( naming_method='convention', export_directory=folder)

def pIMOS_export(rr, archive_dir):

    folder = pIMOS_get_my_export_folder(rr, archive_dir)
    print(folder)

    # rr.folder = folder
    # rr.file_ = serial
    
    rr.export( naming_method='convention', export_directory=folder)

def row_to_attrs(row):
    """
    Alias for autonomous_row_to_attrs
    """

    return autonomous_row_to_attrs(row)

def autonomous_row_to_attrs(row):
    """
    Take a row from the autonomous instrument metadata log and return a pIMOS-ready attributes dict. 
    """
    attributes =  {
        'project': row['Project'],
        'trip': row['Trip Recovered'],
        'site': row['Site'],
        'site_station': row['StationID'],
        'nominal_latitude': float(row['Latitude']),
        'nominal_longitude': float(row['Longitude']),
        'nominal_site_depth': float(row['Depth']),
        'nominal_instrument_height_asb': float(row['InstrumentHeight']),
        'nominal_instrument_orientation': row['Orientation'],
        'timezone': row['TimeZone'],
        'instrument_model': row['InstrumentType'],
        'instrument_serial_number': row['SerialNo'],
        }

    if row['Type'].lower() == 'moored':
        attributes['is_profile_data'] = 0
    elif row['Type'].lower() in['profile', 'profileaux' ]:
        attributes['is_profile_data'] = 1
    else:
        warnings.warn('Unrecognised enntry in Type column: {}'.format(row['Type'].lower()) )

    for att in attributes:
        if isinstance(attributes[att], str):
            attributes[att] = attributes[att].strip()
        elif np.isnan(attributes[att]):
            attributes[att] = ''

    return attributes

    
def profiler_row_to_attrs(row):
    """
    Take a row from the profiling instrument metadata log and return a pIMOS-ready attributes dict. 
    """
    attributes = {
        'project': row['Project'],
        'raw_file_name': row['Filename'],
        'site_station': row['StationID'],
        'trip': row['Trip'],
        'is_profile_data': 1,
        'timezone': 'UTC',
        'instrument_model':  row['InstrumentType'],
        'nominal_latitude': row['Profiler in [Lat]'],
        'nominal_longitude': row['Profiler in [Lon]'],
        'nominal_site_depth': row['Total Depth (m)']}

    for att in attributes:
        if isinstance(attributes[att], str):
            attributes[att] = attributes[att].strip()
        elif np.isnan(attributes[att]):
            attributes[att] = ''

    return attributes


def nonempty_attrs(rr):
    md = {k: v for k, v in rr.ds.attrs.items() if v}
    return {k: md[k] for k in md if not pd.isna(md[k])}