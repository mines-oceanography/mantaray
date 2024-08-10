
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
) -> xr.DataSet:
    """Propagate a single ray without considering the effect of currents

    Parameters
    ----------

    Return
    ------
    xr.Dataset :
        A dataset containing the time evolution of the ray

    Examples
    --------
    >>> mantaray.single_ray(-1000, 0, 0.01, 0, 10, 2, "island.nc")
    """
    tmp = _mantaray.single_ray(
        x0, y0, kx0, ky0, duration, step_size, bathymetry, current
    )

    tmp = np.array(tmp)
    varnames = ["time", "x", "y", "kx", "ky"]
    output = xr.Dataset(
        data_vars = {v:(["time_step"],t) for (v,t) in zip(varnames, tmp.T)},
        attrs = {
            "date_created": str(datetime.datetime.now()),
        }
    )

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
) -> xr.DataSet:
    """Ray tracing for multiple initial conditions

    For a given set of initial conditions, progapage those multiple rays in
    parallel and return the projections for each ray

    Parameters
    ----------
    x0 : Sequence[float]
    y0 : Sequence[float]
    kx0 : Sequence[float]
    ky0 : Sequence[float]
    duration : float
    step_size : float
    bathymetry : str
    current : str

    Returns
    -------
    xr.Dataset :
        A dataset containing the time evolution of multiple rays.
    """
    tmp = _mantaray.ray_tracing(
        x0, y0, kx0, ky0, duration, step_size, bathymetry, current
    )

    tmp = np.array(tmp)

    varnames = ["time", "x", "y", "kx", "ky"]
    output = xr.Dataset(
        data_vars = {v:(["time_step", "ray"],t) for (v,t) in zip(varnames, tmp.T)},
        attrs = {
            "date_created": str(datetime.datetime.now()),
        }
    ).transpose("ray", "time_step")

    return output
