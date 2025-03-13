"""Demo for testing the core module. We shall improve this."""

import mantaray

import numpy as np
import xarray as xr


def deep_water_constant_depth():
    ds = xr.Dataset(
        data_vars=dict(
            depth=(["x", "y"], 10_000 * np.ones((3, 3))),
        ),
        coords=dict(
            x=("x", [-1e4, 0, 1e4]),
            y=("y", [-1e4, 0, 1e4]),
        ),
        attrs=dict(description="Test dataset, "),
    )

    return ds


def zero_current_field():
    ds = xr.Dataset(
        data_vars=dict(
            u=(["x", "y"], 0.01 * np.ones((3, 3))),
            v=(["x", "y"], 0.01 * np.ones((3, 3))),
        ),
        coords=dict(
            x=("x", [-1e8, 0, 1e8]),
            y=("y", [-1e8, 0, 1e8]),
        ),
        attrs=dict(description="Test dataset, zero velocity field"),
    )

    return ds


def test_single_ray(tmp_path):
    ds = deep_water_constant_depth()
    ds.to_netcdf(tmp_path / "island.nc", format="NETCDF3_CLASSIC")

    ds = zero_current_field()
    ds.to_netcdf(tmp_path / "current.nc", format="NETCDF3_CLASSIC")

    ds = mantaray.single_ray(
        -1000, 0, 0.01, 0, 10, 2, tmp_path / "island.nc", tmp_path / "current.nc"
    )

    assert ds.sizes["time_step"] == 6
    assert (ds.kx == 0.01).all()
    assert (ds.ky == 0.0).all()


def test_multiple_rays(tmp_path):
    """Test multiple rays."""
    ds = deep_water_constant_depth()
    ds.to_netcdf(tmp_path / "island.nc", format="NETCDF3_CLASSIC")

    ds = zero_current_field()
    ds.to_netcdf(tmp_path / "current.nc", format="NETCDF3_CLASSIC")

    ds = mantaray.ray_tracing(
        3 * [-1000],
        3 * [0],
        3 * [0.01],
        3 * [0],
        10,
        2,
        str(tmp_path / "island.nc"),
        str(tmp_path / "current.nc"),
    )

    assert ds.sizes["time_step"] == 6
    assert ds.sizes["ray"] == 3
    assert (ds.kx == 0.01).all()
    assert (ds.ky == 0.0).all()


def test_rays_variable_length(tmp_path):
    """Test multiple rays with different sizes

    Assumes that in one direction it will reach the boundary first.
    """
    ds = deep_water_constant_depth()
    ds.to_netcdf(tmp_path / "island.nc", format="NETCDF3_CLASSIC")

    ds = zero_current_field()
    ds.to_netcdf(tmp_path / "current.nc", format="NETCDF3_CLASSIC")

    ds = mantaray.ray_tracing(
        2 * [-1e3],
        2 * [0],
        [-0.01, 0.01],
        2 * [0],
        1e6,
        20,
        str(tmp_path / "island.nc"),
        str(tmp_path / "current.nc"),
    )
