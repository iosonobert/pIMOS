import xarray as _xr, warnings
import imp

import pIMOS.xrwrap.seabird_37_39_56 as seabird_37_39_56
import pIMOS.xrwrap.wetlabs_ntu as wetlabs_ntu

try:
    import pIMOS.xrwrap.nortek_signature as nortek_signature
    import pIMOS.xrwrap.nortek_vector as nortek_vector
    import pIMOS.xrwrap.rdi_adcp as rdi_adcp
    has_dolfyn = True
except AttributeError:
    print('Need to update dolfyn to version 1')
    has_dolfyn = False

import pIMOS.xrwrap.solander_ctd as solander_ctd
import pIMOS.xrwrap.lisst as lisst

try:
    imp.find_module('pyODAS')
    import pIMOS.xrwrap.rsi_vmp as rsi_vmp
    has_rsi_vmp = True
except ImportError:
    print('pyODAS not found, not importing pIMOS.xrwrap.rsi_vmp')
    has_rsi_vmp = False

import pIMOS.xrwrap.solander_ctd as solander_ctd

import pIMOS.utils.quality_control
import pIMOS.utils.plot
import pIMOS.utils.modify

from ._version import __version__

# import pIMOS.xrwrap.pimoswrap as pimoswrap 

classes = [
        seabird_37_39_56.SEABIRD_37_39_56,
        solander_ctd.SOLANDER_CTD,
        wetlabs_ntu.WETLABS_NTU,
        lisst.LISST,
    ]

del seabird_37_39_56
del wetlabs_ntu


del lisst

if has_dolfyn:
    classes.append(nortek_vector.NORTEK_VECTOR)
    classes.append(nortek_signature.NORTEK_SIGNATURE)
    classes.append(rdi_adcp.RDI_ADCP_PD02)
    del nortek_signature
    del nortek_vector
    del rdi_adcp

if has_rsi_vmp:
    classes.append(rsi_vmp.RSI_VMP)
    del rsi_vmp

del solander_ctd

def load_pimos_nc(filename):
    """
    This will load an nc and retuer a wrapped nc. Does so by looping through all wrappers to find the right one based on the title.   
    
    Inputs:
        filename: fullpath to pimos netcdf file
    """
    
    if type(filename) is list:
        ds = _xr.open_mfdataset(filename)   
        print('Opened {}'.format('multtfile dataset')) 
        ds.attrs['outfile_append'] = '' # THis should be cleared for mf datasets

    else:
        ds = _xr.open_dataset(filename)
        print('Opened {}'.format(filename)) 

    return wrap_pymos_ds(ds)
    
def wrap_pymos_ds(ds):
    """
    This will load an nc and retuer a wrapped nc. Does so by looping through all wrappers to find the right one based on the title.   
    
    Inputs:
        ds: pymos xarray dataset
    """
    
    if not 'title' in ds.attrs:
        raise(Exception('File has no title, not a pIMOS file'))
    if not 'source' in ds.attrs:
        raise(Exception('File has no source, not a pIMOS file'))

    title = ds.attrs['title']
    source = ds.attrs['source']
    print('   Title is "{}"'.format(ds.attrs['title'])) 
    print('   Source is "{}"'.format(ds.attrs['source'])) 
    
    if not source == 'pIMOS':
        warnings.warn('Source should be pIMOS. Setting tp pIMOS and warning for now, will become an error in future.')
        ds.attrs['source'] = 'pIMOS'
    
    for c in classes:
        print(c.class_attrs['title'])
        if c.class_attrs['title'] == title:
            print('Dataset appears to be a {}'.format(c))
            return c(ds)
        
    st = ' | '.join(['{}'.format(c) for c in classes])
    
    raise(Exception('File does not match any pIMOS dataset [Options are: {}] [could return generic here instead of error.]'.format(st)))
        