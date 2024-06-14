
from . import _mantaray


def single_ray(x0: float, y0: float, kx0: float, ky0: float, duration: float, step_size: float, filename: str):
    """Propagate a single ray without considering the effect of currents

    Parameters
    ----------

    Return
    ------

    Examples
    --------
    >>> mantaray.single_ray(-1000, 0, 0.01, 0, 10, 2, "island.nc")
    """
    output = _mantaray.single_ray(x0, y0, kx0, ky0, duration, step_size, filename)
    return output
