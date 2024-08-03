
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
):
    """Propagate a single ray without considering the effect of currents

    Parameters
    ----------

    Return
    ------

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
        data_vars = {v:(["time"],t) for (v,t) in zip(varnames, tmp.T)},
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
):
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
    List[List[tuple[float, float, float, float, float]]]:
        List of multiple rays, one for each initial condition.
    """
    tmp = _mantaray.ray_tracing(
        x0, y0, kx0, ky0, duration, step_size, bathymetry, current
    )

    tmp = np.array(tmp)

    varnames = ["time", "x", "y", "kx", "ky"]
    output = xr.Dataset(
        data_vars = {v:(["time", "trajectory"],t) for (v,t) in zip(varnames, tmp.T)},
        attrs = {
            "date_created": str(datetime.datetime.now()),
        }
    ).transpose("trajectory", "time")

    return output
