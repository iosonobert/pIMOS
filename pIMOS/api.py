"""
pIMOS API

Single source of truth for load_pimos_nc and wrap_pymos_ds.

The `classes` list is injected by pIMOS/__init__.py after it has resolved all
optional dependencies.  Nothing in this module imports pIMOS instrument modules
directly, so missing optional dependencies do not affect users who only need
the loader functions.
"""
import xarray as xr
import warnings

# Populated by pIMOS/__init__.py after optional imports are resolved.
classes = []


def load_pimos_nc(filename):
    """
    Load a pIMOS netCDF file and return the appropriate wrapped dataset.

    Parameters
    ----------
    filename : str or list of str
        Path(s) to pIMOS netCDF file(s).
    """
    if isinstance(filename, list):
        ds = xr.open_mfdataset(filename)
        print('Opened multifile dataset')
        ds.attrs['outfile_append'] = ''  # cleared for mf datasets
    else:
        ds = xr.open_dataset(filename)
        print('Opened {}'.format(filename))
    return wrap_pymos_ds(ds)


def wrap_pymos_ds(ds):
    """
    Wrap a raw xarray Dataset in the matching pIMOS instrument class.

    Parameters
    ----------
    ds : xarray.Dataset
        A dataset previously written by pIMOS (must have 'title' and 'source'
        attributes).
    """
    if 'title' not in ds.attrs:
        raise Exception('File has no title attribute — not a pIMOS file')
    if 'source' not in ds.attrs:
        raise Exception('File has no source attribute — not a pIMOS file')

    title = ds.attrs['title']
    source = ds.attrs['source']
    print('   Title  : "{}"'.format(title))
    print('   Source : "{}"'.format(source))

    if source != 'pIMOS':
        warnings.warn(
            'Source should be "pIMOS". Setting to "pIMOS" for now; '
            'this will become an error in a future release.'
        )
        ds.attrs['source'] = 'pIMOS'

    for c in classes:
        if c.class_attrs['title'] == title:
            print('Dataset matches {}'.format(c))
            return c(ds)

    options = ' | '.join(str(c) for c in classes)
    raise Exception(
        'File title "{}" does not match any registered pIMOS class. '
        'Available classes: [{}]'.format(title, options)
    )
        