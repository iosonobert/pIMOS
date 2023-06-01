import datetime
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
import afloat.plot.plotting as zplot
import warnings

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

def plot_echo(rr, db_data, mooring, attributes, variable='echo', width=65, cmap='magma'):
    
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

    add_mooring_dates(db_data, mooring, plt.gca())
    plt.show()

    return fig


def plot_temp(rr, db_data, mooring, attributes, variable='Temperature',\
              plotraw=True, width=65, experiment=None, recovered=None):
    
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

    if plotraw:
        rr.ds[variable].plot(label='Raw')
    rr.get_qaqc_var(variable).plot(label='QAQC')

    plt.xlim(rr.ds.time.values[[0, -1]])
    if plotraw:
        plt.legend()
    
    plt.grid()

    title= ' | '.join([str(attributes[i]) for i in attributes])
    plt.title(title)
    
    if not rr.attrs['is_profile_data']:
        add_mooring_dates(db_data, mooring, plt.gca(), experiment=experiment, recovered=recovered)
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

def pIMOS_export(rr, archive_dir, file_version, model, file_append=''):

    folder = os.path.join(archive_dir, model)
    folder = os.path.join(archive_dir, 'pimos_v'+__version__)
    if not os.path.exists(folder):
        os.mkdir(folder)

    fv = 'FV{:02.0f}'.format(file_version)
    folder = os.path.join(folder, fv)
    if not os.path.exists(folder):
        os.mkdir(folder)

    folder = os.path.join(folder, model)
        
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