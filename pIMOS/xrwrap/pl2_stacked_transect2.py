# -*- coding: utf-8 -*-
"""
Continuous PL2 transect stacker for Triaxus-style profiling data.
"""

import datetime
from collections import OrderedDict

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
from scipy.spatial import QhullError
import xarray as xr

import pIMOS.xrwrap.pimoswrap as pimoswrap

font = {'weight': 'normal', 'size': 15}
matplotlib.rc('font', **font)


class PL2_STACKED_TRANSECT2(pimoswrap.pimoswrap):

    _default_attrs = pimoswrap.pimoswrap._default_attrs.copy()
    _default_attrs['process_level'] = 2

    class_attrs = {
        'title': 'Continuous transect made by stacking profiling data',
        'source': 'pIMOS',
    }

    def __init__(self, ds):
        print('Initialising accessor.')
        self.ds = ds

        self.store_raw_file_attributes(ds)
        self.enforce_these_attrs(self.class_attrs)


def from_rvpctd_continuous(files_or_datasets, stack_variables, znew, **kwargs):
    """
    Build a continuous z-time PL2 transect from one or more RVPCTD-like datasets.

    Parameters
    ----------
    files_or_datasets : sequence
        Sequence of netcdf paths or xarray datasets.
    stack_variables : sequence of str
        Variables to grid onto the z-time product.
    znew : array-like
        Target z grid.

    Keyword Parameters
    ------------------
    dt_minutes : float, default 5
        Regular output timestep in minutes.
    depth_variable : str, default 'Depth'
        Name of the depth variable in the source datasets.
    time_variable : str, default 'time'
        Name of the time variable/dimension in the source datasets.
    latitude_variable : str, default 'external_latitude'
        Name of the latitude variable in the source datasets.
    longitude_variable : str, default 'external_longitude'
        Name of the longitude variable in the source datasets.
    method : str, default 'linear'
        Interpolation backend. Only 'linear' is implemented for now.
    start : optional datetime-like
        Start time of the output grid. Defaults to first sample.
    end : optional datetime-like
        End time of the output grid. Defaults to last sample.
    """

    dt_minutes = kwargs.pop('dt_minutes', 5)
    depth_variable = kwargs.pop('depth_variable', 'Depth')
    time_variable = kwargs.pop('time_variable', 'time')
    latitude_variable = kwargs.pop('latitude_variable', 'external_latitude')
    longitude_variable = kwargs.pop('longitude_variable', 'external_longitude')
    method = kwargs.pop('method', 'linear')
    start = kwargs.pop('start', None)
    end = kwargs.pop('end', None)

    if kwargs:
        raise TypeError('Unexpected keyword arguments: {}'.format(', '.join(kwargs.keys())))

    datasets, source_files = _load_input_datasets(files_or_datasets)
    ds_merged, joined_attrs = _merge_continuous_inputs(
        datasets,
        depth_variable=depth_variable,
        time_variable=time_variable,
        latitude_variable=latitude_variable,
        longitude_variable=longitude_variable,
        stack_variables=stack_variables,
    )

    time_values = pd.to_datetime(ds_merged[time_variable].values)
    if start is None:
        start = time_values.min()
    else:
        start = pd.to_datetime(start)

    if end is None:
        end = time_values.max()
    else:
        end = pd.to_datetime(end)

    time_grid = _build_time_grid(start, end, dt_minutes)
    znew = np.asarray(znew, dtype=float)

    latitude = np.asarray(ds_merged[latitude_variable].values, dtype=float)
    longitude = np.asarray(ds_merged[longitude_variable].values, dtype=float)
    cumulative_distance = _compute_cumulative_distance(time_values, latitude, longitude)

    output = {}
    for stack_variable in stack_variables:
        output[stack_variable] = _grid_time_depth_variable(
            ds_merged,
            stack_variable=stack_variable,
            depth_variable=depth_variable,
            time_variable=time_variable,
            time_grid=time_grid,
            znew=znew,
            method=method,
        )

    output['latitude'] = _interp_to_time_grid(time_values, latitude, time_grid)
    output['longitude'] = _interp_to_time_grid(time_values, longitude, time_grid)
    output['distance'] = _interp_to_time_grid(time_values, cumulative_distance, time_grid)

    ds_stacked = xr.Dataset(attrs={})
    coords = OrderedDict()
    coords.update({'z': znew})
    coords.update({'time': time_grid.values})

    ds_stacked.update({
        'distance': xr.DataArray(
            output['distance'],
            dims={'time': time_grid.values},
            name='distance',
            attrs={
                'long_name': 'cumulative_distance_along_track',
                'units': 'm',
            },
            coords={'time': time_grid.values},
        )
    })
    ds_stacked.update({
        'latitude': xr.DataArray(
            output['latitude'],
            dims={'time': time_grid.values},
            name='latitude',
            attrs={
                'long_name': 'latitude',
                'standard_name': 'latitude',
                'units': 'degrees_north',
            },
            coords={'time': time_grid.values},
        )
    })
    ds_stacked.update({
        'longitude': xr.DataArray(
            output['longitude'],
            dims={'time': time_grid.values},
            name='longitude',
            attrs={
                'long_name': 'longitude',
                'standard_name': 'longitude',
                'units': 'degrees_east',
            },
            coords={'time': time_grid.values},
        )
    })

    for stack_variable in stack_variables:
        ds_stacked.update({
            stack_variable: xr.DataArray(
                output[stack_variable],
                dims=coords,
                name=stack_variable,
                attrs={},
                coords=coords,
            )
        })

    rr = PL2_STACKED_TRANSECT2(ds_stacked)

    attrs_that_must_be_equal = [
        'project',
        'trip',
        'site',
        'instrument_make',
        'instrument_model',
        'instrument_serial_number',
        'site_station',
        'disclaimer',
        'timezone',
        'is_profile_data',
    ]

    for attr in attrs_that_must_be_equal:
        if attr not in joined_attrs:
            continue
        split_attr = joined_attrs[attr]
        if len(split_attr) == 0:
            continue
        if not np.all(np.array(split_attr) == split_attr[0]):
            print(split_attr)
            raise Exception('All files must have the same {}'.format(attr))
        rr.ds.attrs[attr] = split_attr[0]
        joined_attrs.pop(attr)

    rr.ds.attrs['source'] = ';'.join(source_files) if source_files else 'pIMOS'
    for attr in joined_attrs.keys():
        joined_attrs[attr] = [str(item) for item in joined_attrs[attr]]
        rr.ds.attrs[attr] = ';'.join(joined_attrs[attr])

    rr.ds.attrs['stacking_mode'] = 'continuous_triaxus'
    rr.ds.attrs['interpolation_method'] = method
    rr.ds.attrs['interpolation_time_step_minutes'] = dt_minutes
    rr.ds.attrs['created'] = datetime.datetime.now().isoformat()

    return rr


def _load_input_datasets(files_or_datasets):
    datasets = []
    source_files = []

    for item in files_or_datasets:
        if isinstance(item, xr.Dataset):
            datasets.append(item)
            source_files.append(item.attrs.get('raw_file_name', 'in_memory_dataset'))
        else:
            datasets.append(xr.open_dataset(item))
            source_files.append(str(item))

    if len(datasets) == 0:
        raise ValueError('files_or_datasets must contain at least one dataset or path.')

    return datasets, source_files


def _merge_continuous_inputs(datasets, depth_variable, time_variable, latitude_variable, longitude_variable, stack_variables):
    attrs_to_join = [
        'project',
        'trip',
        'trip_deployed',
        'site',
        'site_station',
        'instrument_make',
        'instrument_model',
        'instrument_serial_number',
        'raw_file_name',
        'raw_file_directory',
        'raw_file_attributes',
        'disclaimer',
        'nominal_latitude',
        'nominal_longitude',
        'nominal_site_depth',
        'timezone',
        'is_profile_data',
    ]

    joined_attrs = {attr: [] for attr in attrs_to_join}
    trimmed_datasets = []
    required_variables = [depth_variable, latitude_variable, longitude_variable] + list(stack_variables)

    for ds in datasets:
        if time_variable not in ds.dims:
            raise ValueError('Input dataset must have {} as a dimension.'.format(time_variable))

        for variable in required_variables:
            if variable not in ds:
                raise ValueError('Input dataset missing required variable {}.'.format(variable))

        keep_vars = [time_variable] + required_variables
        keep_vars = [name for name in keep_vars if name in ds]
        trimmed = ds[keep_vars].copy()
        trimmed = trimmed.sortby(time_variable)
        trimmed_datasets.append(trimmed)

        for attr in attrs_to_join:
            if attr in ds.attrs:
                joined_attrs[attr].append(ds.attrs[attr])

    ds_merged = xr.concat(trimmed_datasets, dim=time_variable)
    ds_merged = ds_merged.sortby(time_variable)

    _, unique_index = np.unique(pd.to_datetime(ds_merged[time_variable].values).view('int64'), return_index=True)
    unique_index = np.sort(unique_index)
    ds_merged = ds_merged.isel({time_variable: unique_index})

    return ds_merged, joined_attrs


def _build_time_grid(start, end, dt_minutes):
    freq_seconds = int(round(float(dt_minutes) * 60.0))
    if freq_seconds <= 0:
        raise ValueError('dt_minutes must be positive.')

    if end < start:
        raise ValueError('end must be greater than or equal to start.')

    return pd.date_range(start=start, end=end, freq='{}s'.format(freq_seconds))


def _interp_to_time_grid(native_time, values, time_grid):
    native_time = pd.to_datetime(native_time)
    values = np.asarray(values, dtype=float)

    valid = (~pd.isna(native_time)) & np.isfinite(values)
    native_ns = native_time.view('int64')[valid]
    values = values[valid]

    if len(native_ns) < 2:
        return np.full((len(time_grid),), np.nan)

    sort_idx = np.argsort(native_ns)
    native_ns = native_ns[sort_idx]
    values = values[sort_idx]

    native_ns, unique_index = np.unique(native_ns, return_index=True)
    values = values[unique_index]

    return np.interp(time_grid.view('int64'), native_ns, values, left=np.nan, right=np.nan)


def _compute_cumulative_distance(native_time, latitude, longitude):
    latitude_interp = _interp_to_time_grid(native_time, latitude, native_time)
    longitude_interp = _interp_to_time_grid(native_time, longitude, native_time)

    distance = np.full((len(latitude_interp),), np.nan)
    valid = np.isfinite(latitude_interp) & np.isfinite(longitude_interp)
    if np.sum(valid) < 2:
        return distance

    valid_index = np.where(valid)[0]
    lat_valid = latitude_interp[valid_index]
    lon_valid = longitude_interp[valid_index]

    dlat = np.diff(lat_valid)
    dlon = np.diff(lon_valid)
    lat_ref = np.deg2rad((lat_valid[1:] + lat_valid[:-1]) / 2.0)
    dy = dlat * 111310.0
    dx = dlon * 111310.0 * np.cos(lat_ref)

    cumulative_valid = np.concatenate(([0.0], np.cumsum(np.sqrt(dx**2 + dy**2))))
    distance[valid_index] = cumulative_valid

    if valid_index[0] > 0:
        distance[:valid_index[0]] = np.nan
    if valid_index[-1] < len(distance) - 1:
        distance[valid_index[-1] + 1:] = np.nan

    return distance


def _grid_time_depth_variable(ds, stack_variable, depth_variable, time_variable, time_grid, znew, method):
    if method != 'linear':
        raise NotImplementedError('Interpolation method {} is not implemented yet.'.format(method))

    data = np.asarray(ds[stack_variable].values, dtype=float)
    depth = np.asarray(ds[depth_variable].values, dtype=float)
    time_vals = pd.to_datetime(ds[time_variable].values)

    if ds[stack_variable].dims != (time_variable,):
        raise ValueError('Variable {} must only depend on {}.'.format(stack_variable, time_variable))

    valid = (~pd.isna(time_vals)) & np.isfinite(depth) & np.isfinite(data)
    if np.sum(valid) < 3:
        return np.full((len(znew), len(time_grid)), np.nan)

    time_num = time_vals.view('int64').astype(float)
    target_time_num = time_grid.view('int64').astype(float)

    points = np.column_stack((time_num[valid], depth[valid]))
    grid_time, grid_z = np.meshgrid(target_time_num, znew)

    try:
        return interpolate.griddata(
            points,
            data[valid],
            (grid_time, grid_z),
            method='linear',
            fill_value=np.nan,
        )
    except QhullError:
        # Degenerate point clouds (e.g. near-collinear in time-depth space) can
        # make linear Delaunay triangulation fail. Retry in a scaled space with
        # tiny deterministic jitter; if still failing, fall back to nearest.
        t = points[:, 0]
        z = points[:, 1]

        t_scale = np.nanstd(t)
        z_scale = np.nanstd(z)
        if (not np.isfinite(t_scale)) or (not np.isfinite(z_scale)) or t_scale == 0.0 or z_scale == 0.0:
            return interpolate.griddata(
                points,
                data[valid],
                (grid_time, grid_z),
                method='nearest',
                fill_value=np.nan,
            )

        t0 = np.nanmean(t)
        z0 = np.nanmean(z)
        points_scaled = np.column_stack(((t - t0) / t_scale, (z - z0) / z_scale))

        eps = 1e-10
        idx = np.arange(points_scaled.shape[0], dtype=float)
        jitter = np.column_stack((np.sin(idx), np.cos(idx))) * eps
        points_scaled = points_scaled + jitter

        grid_time_scaled = (grid_time - t0) / t_scale
        grid_z_scaled = (grid_z - z0) / z_scale

        try:
            return interpolate.griddata(
                points_scaled,
                data[valid],
                (grid_time_scaled, grid_z_scaled),
                method='linear',
                fill_value=np.nan,
            )
        except QhullError:
            return interpolate.griddata(
                points,
                data[valid],
                (grid_time, grid_z),
                method='nearest',
                fill_value=np.nan,
            )