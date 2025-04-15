import math
import numpy as np
import cmocean
import matplotlib.pyplot as plt
import matplotlib.animation as animation

g = 9.81 # Acceleration due to gravity [m/s^2]

def period2wavenumber(T):
    """
    Convert wave period to wavenumber for deep water waves.

    Parameters
    ----------
    T : float
        Wave period in seconds.

    Returns
    -------
    k : float
        Wavenumber in radians per meter.
    """
    k = (2 * math.pi)**2 / (g * T**2)
    return k


def group_velocity(k):
    """
    Compute the group velocity for deep water waves given a wavenumber.

    Parameters
    ----------
    k : float
        Wavenumber in radians per meter.

    Returns
    -------
    c_g : float
        Group velocity in meters per second.
    """
    c_g = (g / k)**0.5 / 2
    return c_g


def compute_cfl(x, y, k0):
    """
    Compute the optimal time step for numerical modeling based on CFL condition.

    Parameters
    ----------
    x : np.ndarray
        1D array of x-coordinate values (in meters).
    y : np.ndarray
        1D array of y-coordinate values (in meters).
    k0 : float
        Initial wavenumber in radians per meter.

    Returns
    -------
    cfl : float
        CFL-based time step in seconds.
    """
    dd = np.min([np.diff(x).mean(), np.diff(y).mean()])
    c_g = group_velocity(k0)
    cfl = dd / c_g
    return cfl


def compute_duration(x, k0):
    """
    Estimate model duration based on domain size and wave group velocity.

    Parameters
    ----------
    x : np.ndarray
        1D array of x-coordinate values (in meters).
    k0 : float
        Initial wavenumber in radians per meter.

    Returns
    -------
    duration : int
        Estimated model duration in seconds.
    """
    c_g = group_velocity(k0)
    return round(x.max() / c_g)


def animate_rays(X, Y, background, bando, style, ray_sample=1, time_sample=10):
    """
    Create an animation of ray paths over a background field.

    Parameters
    ----------
    X : np.ndarray
        2D array of x-coordinates for background field.
    Y : np.ndarray
        2D array of y-coordinates for background field.
    background : np.ndarray
        2D array of background values (e.g., speed or depth).
    bando : xr.Dataset
        Dataset containing ray trajectory information with dimensions
        'ray', 'time_step', and variables 'x' and 'y'.
    style : str
        Style of background field: 'currents' or 'bathymetry'.
    ray_sample : int, optional
        Interval for subsampling rays (default is 1).
    time_sample : int, optional
        Interval for subsampling animation time steps (default is 10).

    Returns
    -------
    anim : matplotlib.animation.FuncAnimation
        Animation object that can be saved or displayed.
    """
    time_steps = bando.time_step.size

    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)

    if style == 'currents':
        cf = ax.contourf(X, Y, background, cmap=cmocean.cm.speed, levels=50)
        cbar = fig.colorbar(cf)
        cbar.set_label("Speed [m/s]")
    if style == 'bathymetry':
        cf = ax.contourf(X, Y, background, cmap=cmocean.cm.deep, levels=50)
        cbar = fig.colorbar(cf)
        cbar.set_label("Depth [m]")

    ray_lines = []
    for i in range(0, bando.ray.size, ray_sample):
        color = 'black' if style == 'currents' else 'white'
        ray, = ax.plot([], [], lw=0.78, color=color)
        ray_lines.append(ray)

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title("Ray Tracing Animation")

    def animate(frame):
        for i, ray_line in enumerate(ray_lines):
            ray = bando.isel(ray=i * ray_sample).sel(time_step=slice(0, frame))
            ray_line.set_data(ray.x, ray.y)
        ax.set_title(f"Ray Tracing Animation - Time Step {frame}")

    anim = animation.FuncAnimation(fig, animate, frames=range(0, time_steps, time_sample), interval=100)

    plt.close(fig)
    return anim