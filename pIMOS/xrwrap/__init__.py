"""Public compatibility facade for registered pIMOS wrappers.

Direct submodule imports such as ``import pIMOS.xrwrap.wetlabs_ntu`` continue to
work because the submodules still live in this package. Higher-level attribute
access such as ``pIMOS.xrwrap.wetlabs_ntu`` is routed through ``pIMOS.wrappers``
so unavailable wrappers raise a clear dependency error on access.
"""


def __getattr__(name):
    if name == 'register_wrappers':
        from .._private_wrappers import register_wrappers
        return register_wrappers

    import pIMOS
    return getattr(pIMOS.wrappers, name)


def __dir__():
    import pIMOS
    return sorted(set(globals()) | set(dir(pIMOS.wrappers)))
