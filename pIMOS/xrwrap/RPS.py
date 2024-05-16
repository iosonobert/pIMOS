import os
import pIMOS.xrwrap.xrwrap as xrwrap
from pIMOS.read import RPS as read_RPS
import matplotlib.pyplot as plt


# Main read function
def from_nc(filename, trip=None, RPSvars=None, parse_times=True, verbose=True):
    '''
    Main file reading function for RPS format netCDF files for pIMOS

    Parameters
    ----------
    filename : str
        The file location (full path)
    RPSvars : dict
        Dictionary to specify which RPS variables to capture and 
        how to rename them. Passed through to parse func.
    parse_times : bool, optional
        A flag used to force datetime64 parsing (default is
        True). This flag indicates whether the file is in 
        the new or old RPS format.
    verbose : bool, optional
        A flag to set printing info on or off (default is
        True (on)).

    Returns
    -------
    rr
        pIMOS xwrap class
    ds
        xrray DataSet with pIMOS attributes
    dsr
        xarray DataSet with original attributes (raw load)
    '''

    # Filename checker
    folder, file = xrwrap.parse_infile(filename, verbose=verbose)
    fullpath = os.path.join(folder, file)

    # Call RPS reader
    if RPSvars:
        dsr, ds = read_RPS.process_rps_file(fullpath, RPSvars=RPSvars,\
                                            parse_times=parse_times, verbose=verbose)
    else:
        dsr, ds = read_RPS.process_rps_file(fullpath,\
                                            parse_times=parse_times, verbose=verbose)
                
    # xwrap RPS class
    rr = RPS(ds, verbose=verbose) # change back later

    # add key attrs
    rr.ds.attrs['raw_file_name']      = os.path.split(filename)[1]
    rr.ds.attrs['raw_file_directory'] = os.path.split(filename)[0]
    rr.ds.attrs['site_station']       = read_RPS.lscm(dsr.attrs['location'])
    rr.ds.attrs['instrument_model']   = read_RPS.lscm(dsr.attrs['logging_system'])

    if len(rr.ds.attrs['instrument_serial_number']) == 0:
        try:
            rr.ds.attrs['instrument_serial_number'] = dsr.attrs['dms_series_id']
            if verbose:
                print('Setting serial number to RPS dataset ID')
        except:
            if verbose:
                print('No serial number')
            pass

    # Add trip (important for multi-phase deployments)
    if trip is not None:
        rr.ds.attrs['trip'] = trip
    else:
        # Add a trip as the month file starts
        rr.ds.attrs['trip'] = read_RPS.lscm(str(rr.ds['time'][0].values.astype('datetime64[M]')))

    if verbose:
        print('Done')
    return rr, ds, dsr




class RPS(xrwrap.xrwrap):
    
    class_attrs = {
            'title': 'Measured data parsed from an RPS data file.',
            'source': 'pIMOS RPS' # Could be more specific.
        }

    def __init__(self, ds, verbose=True):
        # Not sure when to invole verbose
        self.verbose = verbose

        if self.verbose:
            print('Initialising accessor.')
        self.ds = ds # XRWRAP compatibility

        self.store_raw_file_attributes(ds)
        
        self.enforce_these_attrs(self.class_attrs)


    def plot_RPS_export(self, plt_vars, labels='short'):

        fig, ax = plt.subplots(len(plt_vars),1, figsize=(14,1.0*len(plt_vars)), gridspec_kw={'hspace':0.07})
        if not hasattr(ax, '__len__'):
            ax = [ax]
        for x, dvar in zip(ax, plt_vars):
            self.ds[dvar].plot(ax=x, lw=0.8)
            x.set_xlim(self.ds['time'][0], self.ds['time'][-1])
            x.grid()
            if labels=='long':
                if 'long_name' in self.ds[dvar].attrs:
                    x.set_ylabel(self.ds[dvar].long_name + '\n[' + self.ds[dvar].units + ']')
                else:
                    x.set_ylabel(dvar + '\n[' + self.ds[dvar].units + ']')
            else:
                x.set_ylabel(dvar)
            if x != ax[-1]:
                x.set_xticklabels('')
            else:
                x.set_xticks(x.get_xticks()) # this avoids a mpl warning
                x.set_xticklabels(x.get_xticklabels(), rotation=0, ha="center")
            if x != ax[0]:
                x.set_title('')
            x.set_xlabel('')
        return fig