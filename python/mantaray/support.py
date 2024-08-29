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
    bathymetry: Optional[str] = None,
    current: Optional[str] = None,
    spacing: int = 5,
) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]:
    '''Plot the ray bundle using matplotlib. Optionally plot bathymetry or current.

    Parameters
    ----------
    ray_bundle : xarray.Dataset
        An xarray dataset containing the rays, time steps, x, y, kx, and ky,
        values

    bathymetry : str, optional
        The path to the bathymetry file to plot

    current : str, optional
        The path of the current file to plot

    spacing : int, default=5
        The spacing between plotting nearby current vectors. The default is 5,
        which means every fifth current vector will be plotted.

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
        plot_current(current, ax, spacing)

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
    ds = xr.open_dataset(bathymetry)
    x = ds.x.values
    y = ds.y.values
    h = ds.depth.values
    im = ax.imshow(h, extent=(x[0], x[-1], y[0], y[-1]))
    # Color bar
    colorbar = ax.figure.colorbar(im, ax=ax)
    colorbar.ax.set_ylabel("Depth [m]", rotation=-90, va="bottom")


def plot_bathymetry_contour(bathymetry, ax):
    ds = xr.open_dataset(bathymetry)
    x = ds.x.values
    y = ds.y.values
    h = ds.depth.values
    cs = ax.contour(x, y, h, extent=(x[0], x[-1], y[0], y[-1]), colors="k")
    ax.clabel(cs, cs.levels, inline=True, fontsize=10, colors="k")


def plot_current(current, ax, spacing):
    ds = xr.open_dataset(current)
    x = ds.x.values
    y = ds.y.values
    u = ds.u.values
    v = ds.v.values

    speed = np.sqrt(u**2 + v**2)
    im = ax.imshow(speed, extent=(x[0], x[-1], y[0], y[-1]))
    colorbar = ax.figure.colorbar(im, ax=ax)
    colorbar.ax.set_ylabel("Current Speed [m/s]", rotation=-90, va="bottom")
    X, Y = np.meshgrid(x, y)

    ax.quiver(
        X[::spacing, ::spacing],
        Y[::spacing, ::spacing],
        u[::spacing, ::spacing],
        v[::spacing, ::spacing],
        scale=10,
    )
