from typing import Optional
import matplotlib.axes
import matplotlib.figure
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr


def plot_ray_tracing(
    ray_bundle: xr.Dataset,
    *,
    bathymetry: Optional[xr.Dataset] = None,
    current: Optional[xr.Dataset] = None,
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    '''Plot the ray bundle using matplotlib. Optionally plot bathymetry or current.

    Parameters
    ----------
    ray_bundle : xarray.Dataset
        An xarray dataset containing the rays, time steps, x, y, kx, and ky,
        values

    bathymetry : xarray.Dataset, optional
        The path to the bathymetry file to plot

    current : xarray.Dataset, optional
        The path of the current file to plot

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        The created figure and axes. See matplotlib documentation for
        `plt.subplots()` return type.

    '''
    fig, ax = plt.subplots()

    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Plot of Ray Tracing")

    if current:
        plot_current(current, ax)

    if bathymetry:
        # plot as contour or heatmap?
        if current:
            pass
            # plot_bathymetry_contour(bathymetry, ax)
        else:
            plot_bathymetry_heatmap(bathymetry, ax)

    if ray_bundle:
        plot_ray_bundle(ray_bundle, ax)

    plt.show()
    return fig, ax


def plot_ray_bundle(ray_path, ax):
    number_of_rays = ray_path.ray.size
    for i in range(number_of_rays):
        x = ray_path.x.values[i]
        y = ray_path.y.values[i]
    ax.plot(x, y, color="r")


def plot_bathymetry_heatmap(bathymetry, ax):
    #FIXME: make xarray style
    x = bathymetry.x.values
    y = bathymetry.y.values
    h = bathymetry.depth.values
    im = ax.imshow(h, extent=(x[0], x[-1], y[0], y[-1]))
    # Color bar
    colorbar = ax.figure.colorbar(im, ax=ax)
    colorbar.ax.set_ylabel("Depth [m]", rotation=-90, va="bottom")


def plot_bathymetry_contour(bathymetry, ax):
    #FIXME: make xarray style
    x = bathymetry.x.values
    y = bathymetry.y.values
    h = bathymetry.depth.values
    cs = ax.contour(x, y, h, extent=(x[0], x[-1], y[0], y[-1]), colors="k")
    ax.clabel(cs, cs.levels, inline=True, fontsize=10, colors="k")


def plot_current(current, ax):
    #FIXME: make xarray style
    x = current.x.values
    y = current.y.values
    u = current.u.values
    v = current.v.values

    speed = np.sqrt(u**2 + v**2)
    im = ax.imshow(speed, extent=(x[0], x[-1], y[0], y[-1]))
    colorbar = ax.figure.colorbar(im, ax=ax)
    colorbar.ax.set_ylabel("Current Speed [m/s]", rotation=-90, va="bottom")
    X, Y = np.meshgrid(x, y)

    ax.quiver(
        X,
        Y,
        u,
        v,
        scale=10,
    )
