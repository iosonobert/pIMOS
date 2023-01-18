import xarray as _xr, warnings

import pIMOS.xrwrap.seabird_37_39_56 as seabird_37_39_56
import pIMOS.xrwrap.wetlabs_ntu as wetlabs_ntu
import pIMOS.xrwrap.nortek_signature as nortek_signature
import pIMOS.xrwrap.nortek_vector as nortek_vector
import pIMOS.xrwrap.lisst as lisst
import pIMOS.xrwrap.rdi_adcp as rdi_adcp
import pIMOS.xrwrap.rsi_vmp as rsi_vmp
import pIMOS.xrwrap.solander_ctd as solander_ctd

import pIMOS.utils.quality_control
import pIMOS.utils.plot
import pIMOS.utils.modify


# import pIMOS.xrwrap.pimoswrap as pimoswrap 

classes = [
        seabird_37_39_56.SEABIRD_37_39_56,
        rsi_vmp.RSI_VMP,
        solander_ctd.SOLANDER_CTD,
        rdi_adcp.RDI_ADCP_PD02,
        wetlabs_ntu.WETLABS_NTU,
        lisst.LISST,
        nortek_vector.NORTEK_VECTOR,
        nortek_signature.NORTEK_SIGNATURE
    ]

del seabird_37_39_56
del wetlabs_ntu
del nortek_signature
del nortek_vector
del lisst
del rdi_adcp
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
        