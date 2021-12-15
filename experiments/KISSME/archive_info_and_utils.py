import datetime
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
import zutils.plotting as zplot

"""
Various tools that should be pIMOS tools, but that I have not officially found a spot for yet. 
"""

# """
# I have basically made these up for now. 
# """

# moorings = {
#     'SP250': {
#         'start': datetime.datetime(2017, 4, 1, 9, 0, 0),
#         'end': datetime.datetime(2017, 5, 21, 23, 30, 0)
#     },
#     'SP250L': {
#         'start': datetime.datetime(2017, 4, 1, 9, 0, 0),
#         'end': datetime.datetime(2017, 5, 21, 23, 30, 0)
#     },
#     'WP250': {
#         'start': datetime.datetime(2017, 4, 2, 9, 0, 0),
#         'end': datetime.datetime(2017, 5, 21, 23, 30, 0)
#     },
#     'NP250': {
#         'start': datetime.datetime(2017, 4, 2, 9, 0, 0),
#         'end': datetime.datetime(2017, 5, 21, 23, 30, 0)
#     },
#     'FAKE': {
#         'start': datetime.datetime(2017, 4, 2, 9, 0, 0),
#         'end': datetime.datetime(2017, 5, 21, 23, 30, 0)
#     },
# }

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
    if 'deployment_metadata' in db_config:
        db_data['deployment_metadata'] = parse_deployment_metadata(db_config)
    else:
        raise(Exception('Config file incomplete'))
        
    if 'possible_mooring_dates' in db_config:
        db_data['possible_mooring_dates'] = parse_possible_mooring_dates(db_config, recovered=None)
    else:
        raise(Exception('Config file incomplete'))
        
    return db_data

def parse_deployment_metadata(db_config):

    file = os.path.join(db_config['db_root'], db_config['deployment_metadata'])
    df = pd.read_csv(file)
    
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

def plot_temp(rr, db_data, mooring, title, variable='Temperature', plotraw=True, width=65):
    
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
    plt.title(title)
    
    add_mooring_dates(db_data, mooring, plt.gca())
    plt.show()
    
    return fig

def add_mooring_dates(db_data, mooring, ax):
    
    yl = ax.get_ylim()
    
    df = db_data['possible_mooring_dates']

    # Hard Code This For Now
    df = df.loc[strcmpi(df['Recovered'].values, 'kissme_recovery')]
    df = df.loc[strcmpi(df['Experiment'].values, 'kissme')]
    
    df = df.loc[strcmpi(df['Mooring'].values, mooring)]

    ax.plot([df['StartDate'].values[0]]*2, yl, 'r--')
    ax.plot([df['EndDate'].values[0]]*2, yl, 'r--')
    
    pass

def pIMOS_export(rr, archive_dir, model, serial, csv=True):
    
    folder = os.path.join(archive_dir, model)
    if not os.path.exists(folder):
        os.mkdir(folder)

    print(folder)

    rr.folder = folder
    rr.file_ = serial
    if csv:
        rr.export()
    else:
        rr.export(csv=csv)
        