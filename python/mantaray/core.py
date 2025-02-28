import datetime

import numpy as np
import xarray as xr

from . import _mantaray


def single_ray(
    x0: float,
    y0: float,
    kx0: float,
    ky0: float,
    duration: float,
    step_size: float,
    bathymetry: str,
    current: str,
) -> xr.Dataset:
    """Propagate a single ray without considering the effect of currents

    Parameters
    ----------
    x0 : float
        Initial x position of the ray
    y0 : float
        Initial y position of the ray
    kx0 : float
        Initial wavenumber, x component
    ky0 : float
        Initial wavenumber, y component
    duration : float
        Duration of the simulation
    step_size : float
        Time step for the simulation
    bathymetry : str
        Path to a netCDF file containing the bathymetry file. It is expected
        to have x and y dimensions as floats and depth (x, y) as a float,
        where depth is zero at surface and positive downwards.
    current : str
        Paths to a netCDF file containing the current field. It is expected
        to have x and y dimensions as floats and u(x,y) and v(x,y) as floats.

    Return
    ------
    xr.Dataset :
        A dataset containing the time evolution of the ray

    Examples
    --------
    >>> mantaray.single_ray(-1000, 0, 0.01, 0, 10, 2, "island.nc")
    """
    tmp = _mantaray.single_ray(
        x0, y0, kx0, ky0, duration, step_size, str(bathymetry), str(current)
    )

    tmp = np.array(tmp)
    varnames = ["time", "x", "y", "kx", "ky"]
    output = xr.Dataset(
        data_vars={v: (["time_step"], t) for (v, t) in zip(varnames, tmp.T)},
        attrs={
            "date_created": str(datetime.datetime.now()),
        },
    ).set_coords(["time", "x", "y"])

    return output


def ray_tracing(
    x0,
    y0,
    kx0,
    ky0,
    duration: float,
    step_size: float,
    bathymetry: str,
    current: str,
) -> xr.Dataset:
    """Ray tracing for multiple initial conditions

    For a given set of initial conditions, progapage those multiple rays in
    parallel and return the projections for each ray

    Parameters
    ----------
    x0 : Sequence[float]
        Initial x position of the ray
    y0 : Sequence[float]
        Initial y position of the ray
    kx0 : Sequence[float]
        Initial wavenumber, x component
    ky0 : Sequence[float]
        Initial wavenumber, y component
    duration : float
        Duration of the simulation
    step_size : float
        Time step for the simulation
    bathymetry : str
        Path to a netCDF file containing the bathymetry file. It is expected
        x and y dimensions as floats and depth (x, y) as a float, where
        depth is zero at surface and positive downwards.
    current : str
        Paths to a netCDF file containing the current field. It is expected
        to have x and y dimensions as floats and u(x,y) and v(x,y) as floats.

    Returns
    -------
    xr.Dataset :
        A dataset containing the time evolution of multiple rays.
    """
    tmp = _mantaray.ray_tracing(
        x0, y0, kx0, ky0, duration, step_size, bathymetry, current
    )

    varnames = ["time", "x", "y", "kx", "ky"]
    longest_ray = max([len(ray) for ray in tmp])
    bundle = [
        np.pad(ray, ((0, longest_ray - len(ray)), (0, 0)), constant_values=np.nan)
        for ray in tmp
    ]
    ds = xr.Dataset(
        data_vars={
            name: (["time_step", "ray"], v)
            for (name, v) in zip(varnames, np.array(bundle).T)
        }
    )

    ds = ds.set_coords(["time", "x", "y"])
    ds["time_step"] = range(ds.sizes["time_step"])
    ds["ray"] = np.arange(len(ds.ray))
    ds.attrs["date_created"] = str(datetime.datetime.now())

    return ds
