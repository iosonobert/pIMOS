print('INITIALISING pIMOS')
from ._version import __version__
from ._private_wrappers import register_wrappers

wrappers = register_wrappers()
classes = list(wrappers.classes)

_dolfyn_wrapper_names = ('nortek_signature', 'nortek_vector', 'rdi_adcp')
has_dolfyn = all(name in wrappers.available for name in _dolfyn_wrapper_names)
has_rsi_vmp = 'rsi_vmp' in wrappers.available

import pIMOS.utils.quality_control
import pIMOS.utils.plot
import pIMOS.utils.modify

from . import api as _api
_api.classes = classes
from .api import load_pimos_nc, wrap_pymos_ds

__all__ = [
    '__version__',
    'classes',
    'has_dolfyn',
    'has_rsi_vmp',
    'load_pimos_nc',
    'register_wrappers',
    'wrap_pymos_ds',
    'wrappers',
]
