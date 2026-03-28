import numpy as np
import pytest


def _make_minimal_netcdf(path):
    xr = pytest.importorskip("xarray")
    ds = xr.Dataset(
        data_vars={
            "c1": ("time", np.array([1000.0, 1001.0], dtype=float)),
        },
        coords={
            "time": np.array(["2020-01-01T00:00:00", "2020-01-01T00:01:00"], dtype="datetime64[ns]"),
        },
    )
    ds.to_netcdf(path)
    ds.close()


def test_wetlabs_classmethod_from_netcdf(tmp_path):
    wetlabs_mod = pytest.importorskip("pIMOS.xrwrap.wetlabs_ntu")
    WETLABS_NTU = wetlabs_mod.WETLABS_NTU

    nc_path = tmp_path / "wetlabs_test.nc"
    _make_minimal_netcdf(nc_path)

    rr, ds = WETLABS_NTU.from_netcdf(str(nc_path))

    assert isinstance(rr, WETLABS_NTU)
    assert rr.ds is ds
    assert ds.attrs["last_load_file_name"] == nc_path.name
    assert ds.attrs["last_load_directory"] == str(tmp_path)

    ds.close()


def test_wetlabs_module_loader_routes_to_classmethod(tmp_path):
    wetlabs_mod = pytest.importorskip("pIMOS.xrwrap.wetlabs_ntu")
    WETLABS_NTU = wetlabs_mod.WETLABS_NTU
    wetlabs_from_netcdf = wetlabs_mod.from_netcdf

    nc_path = tmp_path / "wetlabs_test_module_entrypoint.nc"
    _make_minimal_netcdf(nc_path)

    rr, ds = wetlabs_from_netcdf(str(nc_path))

    assert isinstance(rr, WETLABS_NTU)
    assert ds.attrs["last_load_file_name"] == nc_path.name

    ds.close()
