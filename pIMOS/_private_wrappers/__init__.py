from dataclasses import dataclass
from importlib import import_module
from types import MappingProxyType


@dataclass(frozen=True)
class WrapperSpec:
    name: str
    module_path: str
    class_name: str | None = None


class UnavailableWrapper:
    def __init__(self, name, module_path, error):
        self.name = name
        self.module_path = module_path
        self.error = error

    @property
    def message(self):
        return str(self.error)

    def as_import_error(self):
        return ImportError(
            "Wrapper '{}' is unavailable from '{}'. {}".format(
                self.name,
                self.module_path,
                self.message,
            )
        )

    def __call__(self, *args, **kwargs):
        raise self.as_import_error() from self.error

    def __getattr__(self, _name):
        raise self.as_import_error() from self.error

    def __repr__(self):
        return "<UnavailableWrapper {}: {}>".format(self.name, self.message)


class WrapperRegistry:
    def __init__(self):
        self._available_modules = {}
        self._available_classes = {}
        self._failed_modules = {}

    def add_available(self, spec, module, wrapper_class=None):
        self._available_modules[spec.name] = module
        if wrapper_class is not None:
            self._available_classes[spec.name] = wrapper_class

    def add_failed(self, spec, error):
        self._failed_modules[spec.name] = UnavailableWrapper(
            spec.name,
            spec.module_path,
            error,
        )

    @property
    def available(self):
        return MappingProxyType(self._available_modules)

    @property
    def failed(self):
        return MappingProxyType(self._failed_modules)

    @property
    def classes(self):
        return tuple(self._available_classes.values())

    @property
    def names(self):
        return tuple(sorted(set(self._available_modules) | set(self._failed_modules)))

    def __contains__(self, name):
        return name in self._available_modules or name in self._failed_modules

    def __getitem__(self, name):
        if name in self._available_modules:
            return self._available_modules[name]
        if name in self._failed_modules:
            raise self._failed_modules[name].as_import_error() from self._failed_modules[name].error
        raise KeyError(name)

    def __getattr__(self, name):
        if name in self._available_modules:
            return self._available_modules[name]
        if name in self._failed_modules:
            raise self._failed_modules[name].as_import_error() from self._failed_modules[name].error
        raise AttributeError(name)

    def __dir__(self):
        return sorted(set(super().__dir__()) | set(self.names))

    def summary(self):
        return {
            'available': sorted(self._available_modules),
            'failed': {name: wrapper.message for name, wrapper in self._failed_modules.items()},
        }


WRAPPER_SPECS = (
    WrapperSpec('rps', 'pIMOS._private_wrappers.RPS', 'RPS'),
    WrapperSpec('drifter', 'pIMOS._private_wrappers.drifter', 'DRIFTER'),
    WrapperSpec('lisst', 'pIMOS._private_wrappers.lisst', 'LISST'),
    WrapperSpec('nortek_signature', 'pIMOS._private_wrappers.nortek_signature', 'NORTEK_SIGNATURE'),
    WrapperSpec('nortek_vector', 'pIMOS._private_wrappers.nortek_vector', 'NORTEK_VECTOR'),
    WrapperSpec('pl2_stacked_mooring', 'pIMOS._private_wrappers.pl2_stacked_mooring', 'PL2_STACKED_MOORING'),
    WrapperSpec('pl2_stacked_transect_cont', 'pIMOS._private_wrappers.pl2_stacked_transect_cont', 'PL2_STACKED_TRANSECT_CONT'),
    WrapperSpec('pl2_stacked_transect_discrete', 'pIMOS._private_wrappers.pl2_stacked_transect_discrete', 'PL2_STACKED_TRANSECT_DISCRETE'),
    WrapperSpec('pl3_stacked_mooring', 'pIMOS._private_wrappers.pl3_stacked_mooring', 'PL3_STACKED_MOORING'),
    WrapperSpec('rbr_duet', 'pIMOS._private_wrappers.rbr_duet', 'RBR_DUET'),
    WrapperSpec('rdi_adcp', 'pIMOS._private_wrappers.rdi_adcp', 'RDI_ADCP_PD02'),
    WrapperSpec('rsi_vmp', 'pIMOS._private_wrappers.rsi_vmp', 'RSI_VMP'),
    WrapperSpec('rvlctd', 'pIMOS._private_wrappers.rvlctd', 'RVLCTD'),
    WrapperSpec('rvpctd', 'pIMOS._private_wrappers.rvpctd', 'RVPCTD'),
    WrapperSpec('seabird_37_39_56', 'pIMOS._private_wrappers.seabird_37_39_56', 'SEABIRD_37_39_56'),
    WrapperSpec('wetlabs_ntu', 'pIMOS._private_wrappers.wetlabs_ntu', 'WETLABS_NTU'),
)


def register_wrappers(specs=WRAPPER_SPECS):
    registry = WrapperRegistry()

    for spec in specs:
        try:
            module = import_module(spec.module_path)
            wrapper_class = getattr(module, spec.class_name) if spec.class_name else None

            # Provide a consistent module-level constructor for blank wrappers.
            # This lets users call pIMOS.wrappers.<name>.blank(...) directly.
            if (
                wrapper_class is not None
                and hasattr(wrapper_class, 'blank')
                and not hasattr(module, 'blank')
            ):
                module.blank = wrapper_class.blank

            registry.add_available(spec, module, wrapper_class)
        except (ImportError, AttributeError) as error:
            registry.add_failed(spec, error)
            # Surface optional dependency problems consistently at import time.
            print(
                "Wrapper '{}' unavailable: {}".format(
                    spec.name,
                    str(error).strip().splitlines()[0],
                )
            )

    return registry
