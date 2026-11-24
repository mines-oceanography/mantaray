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


def test_ray_tracing_arrays_matches_files(tmp_path):
    """Array input gives the same rays as the NetCDF input, and refracts."""
    x = np.arange(0, 20_000.1, 100.0)
    y = np.arange(0, 10_000.1, 100.0)
    depth = np.tile((100 - 90 * y / 10_000)[:, None], (1, len(x)))  # (y, x)
    xr.Dataset(
        dict(depth=(["y", "x"], depth)), coords=dict(x=x, y=y)
    ).to_netcdf(tmp_path / "bathy.nc", format="NETCDF3_CLASSIC")
    zero = np.zeros((3, 3))
    xr.Dataset(
        dict(u=(["y", "x"], zero), v=(["y", "x"], zero)),
        coords=dict(x=[-1e9, 0, 1e9], y=[-1e9, 0, 1e9]),
    ).to_netcdf(tmp_path / "current.nc", format="NETCDF3_CLASSIC")

    args = ([5000.0, 5000.0], [0.0, 0.0], [0.02, 0.03], [0.04, 0.03], 2000.0, 2.0)
    a = mantaray.ray_tracing_arrays(*args, x, y, depth)
    f = mantaray.ray_tracing(
        *args, str(tmp_path / "bathy.nc"), str(tmp_path / "current.nc")
    )
    np.testing.assert_allclose(a.x, f.x, rtol=1e-6, equal_nan=True)
    np.testing.assert_allclose(a.kx, f.kx, rtol=1e-6, equal_nan=True)
    # shoaling: kx is conserved while k grows, so the ray turns towards the normal
    first, last = a.isel(ray=0).dropna("time_step").isel(time_step=[0, -1]).ky.values
    assert last > first
