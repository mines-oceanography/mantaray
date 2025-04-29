import math
import numpy as np
import xarray as xr
import cmocean
import matplotlib.pyplot as plt
import matplotlib.animation as animation

g = 9.81  # Acceleration due to gravity [m/s^2]


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
    k = (2 * math.pi) ** 2 / (g * T**2)
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
    c_g = (g / k) ** 0.5 / 2
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


def animate_rays(X, Y, background, ray_bundle, style, ray_sample=1, time_sample=10):
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
    ray_bundle : xr.Dataset
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
    time_steps = ray_bundle.time_step.size

    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)

    if style == "currents":
        cf = ax.contourf(X, Y, background, cmap=cmocean.cm.speed, levels=50)
        cbar = fig.colorbar(cf)
        cbar.set_label("Speed [m/s]")
    if style == "bathymetry":
        cf = ax.contourf(X, Y, background, cmap=cmocean.cm.deep, levels=50)
        cbar = fig.colorbar(cf)
        cbar.set_label("Depth [m]")

    ray_lines = []
    for i in range(0, ray_bundle.ray.size, ray_sample):
        color = "black" if style == "currents" else "white"
        (ray,) = ax.plot([], [], lw=0.78, color=color)
        ray_lines.append(ray)

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title("Ray Tracing Animation")

    def animate(frame):
        for i, ray_line in enumerate(ray_lines):
            ray = ray_bundle.isel(ray=i * ray_sample).sel(time_step=slice(0, frame))
            ray_line.set_data(ray.x, ray.y)
        ax.set_title(f"Ray Tracing Animation - Time Step {frame}")

    anim = animation.FuncAnimation(
        fig, animate, frames=range(0, time_steps, time_sample), interval=100
    )

    plt.close(fig)
    return anim


def plot_current_field(x_grid, y_grid, ds, skip=75, q_ref=0.5, q_scale=0.1):
    """
    Plot the magnitude and direction of a velocity field using contour and quiver plots.

    Parameters
    ----------
    x_grid : ndarray
        2D array of x-coordinates (meters).
    y_grid : ndarray
        2D array of y-coordinates (meters).
    ds : xarray.Dataset
        Dataset containing 2D velocity components `u` (zonal) and `v` (meridional),
        with the same shape as `x_grid` and `y_grid`.
    skip : int, optional
        Step size for downsampling vectors in the quiver plot to reduce clutter.
        Default is 75 (i.e., plot every 75th vector in each direction).
    q_ref : float, optional
        Reference vector magnitude to display in the quiver key (in m/s).
        Default is 0.5 m/s.
    q_scale : float, optional
        Scale factor for the quiver arrows. Smaller values produce longer arrows.
        Default is 0.1.

    Returns
    -------
    None
        The function displays the plot and returns nothing.

    Notes
    -----
    - Uses `cmocean.cm.speed` for the color map representing velocity magnitude.
    - Quiver vectors are plotted in black and scaled by `q_scale`.
    - Axes are labeled in kilometers, and the plot enforces equal aspect ratio.
    """
    speed = np.sqrt(ds.u**2 + ds.v**2)

    fig, ax = plt.subplots(figsize=(12, 6))

    # Convert x and y to km
    x_km = x_grid / 1000
    y_km = y_grid / 1000

    c = ax.contourf(x_km, y_km, speed, cmap=cmocean.cm.speed, levels=50)
    fig.colorbar(c, ax=ax, label="Speed [m/s]")

    q = ax.quiver(
        x_km[::skip, ::skip],
        y_km[::skip, ::skip],
        ds.u[::skip, ::skip],
        ds.v[::skip, ::skip],
        color="black",
        scale=1 / q_scale,
        width=0.0025,
    )
    ax.quiverkey(q, X=0.9, Y=-0.1, U=q_ref, label=f"{q_ref} [m/s]", labelpos="E")

    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_title("Current Velocity Magnitude and Direction")

    ax.set_aspect("equal")
    ax.grid(linestyle="--")

    plt.show()
    return


def cart2polar(x, y):
    """
    Convert Cartesian coordinates to polar coordinates.

    Parameters
    ----------
    x : array_like
        x-coordinates (meters).
    y : array_like
        y-coordinates (meters).

    Returns
    -------
    r : ndarray
        Radial distance from origin (meters).
    theta : ndarray
        Angle from x-axis in radians (counter-clockwise).
    """
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    return r, theta


def polar2cart_vel(U_theta, theta):
    """
    Convert azimuthal velocity in polar coordinates to Cartesian velocity components.

    Parameters
    ----------
    U_theta : array_like
        Azimuthal (tangential) velocity in polar coordinates (m/s).
    theta : array_like
        Angular coordinate corresponding to each velocity component (radians).

    Returns
    -------
    u : ndarray
        Zonal (x-direction) velocity component (m/s).
    v : ndarray
        Meridional (y-direction) velocity component (m/s).

    Notes
    -----
    Assumes purely azimuthal flow (no radial component).
    """
    u = -U_theta * np.sin(theta)
    v = U_theta * np.cos(theta)
    return u, v


def parabolic_ring_profile(r, r_core, r_outer, U_max=1.0):
    """
    Azimuthal velocity profile for an idealized ring.

    Parameters
    ----------
    r : ndarray
        Radial distance from eddy center (meters).
    r_core : float
        Radius of zero-velocity core (meters).
    r_outer : float
        Outer radius of the eddy (meters).
    U_max : float, optional
        Desired maximum azimuthal velocity (m/s). Default is 1.0 m/s.

    Returns
    -------
    U_theta : ndarray
        Azimuthal velocity at each radial location.
    """
    r0 = 0.5 * (r_core + r_outer)
    width = 0.5 * (r_outer - r_core)

    U_theta = np.zeros_like(r)
    mask = (r >= r0 - width) & (r <= r0 + width)
    rt = r[mask] - r0

    U_norm = (rt + width) * (rt - width)
    U_theta[mask] = -U_max * (U_norm / width**2)

    return U_theta


def generate_parabolic_ring_eddy(
    L_eddy=320_000, U_max=1.0, xv=None, yv=None, core_ratio=0.25
):
    """
    Generate a 2D velocity field representing a parabolic ring eddy.

    Parameters
    ----------
    L_eddy : float, optional
        Size of the square domain in meters (i.e., domain is L_eddy x L_eddy). Default is 320,000 m.
    U_max : float, optional
        Maximum azimuthal velocity in meters per second. Default is 1.0 m/s.
    xv : ndarray, required
        2D array of x-coordinates (meters). If None, a square grid will be generated.
    yv : ndarray, required
        2D array of y-coordinates (meters). Must be the same shape as xv.
    core_ratio : float, optional
        Fraction of the eddy radius used to define the radius of the zero-velocity core.
        For example, core_ratio=0.25 with a 160 km eddy gives a 40 km core. Default is 0.25.

    Returns
    -------
    u : ndarray
        2D array of zonal (east-west) velocity components (m/s).
    v : ndarray
        2D array of meridional (north-south) velocity components (m/s).

    Notes
    -----
    This function models an idealized parabolic ring eddy with a circular velocity profile:
    - The velocity is purely azimuthal (no radial component).
    - The eddy has a central core of zero velocity and a surrounding ring where
      velocity follows a parabolic profile.
    - If xv and yv are not provided, the function will throw an error.
    """
    R_eddy = L_eddy / 2

    r_outer = R_eddy
    r_core = core_ratio * r_outer

    r, theta = cart2polar(xv, yv)

    U_theta = parabolic_ring_profile(r, r_core=r_core, r_outer=r_outer, U_max=U_max)
    u, v = polar2cart_vel(U_theta, theta)

    return u, v

def generate_zonal_jet(U_max, x, y, width):
    """
    Generate a 2D zonal jet velocity field on a regular grid.

    Parameters
    ----------
    U_max : float
        Maximum velocity of the zonal jet (in m/s).
    x : array-like
        1D array of x-coordinates (in meters).
    y : array-like
        1D array of y-coordinates (in meters).
    width : float
        Width (standard deviation) of the Gaussian jet profile (in normalized units).

    Returns
    -------
    ds_jet : xarray.Dataset
        Dataset containing zonal (`u`) and meridional (`v`) velocity components
        defined on a 2D grid, along with coordinate variables `x` and `y`.
    """
    # Normalize y-range to create Gaussian profile across y-direction
    x_profile = np.linspace(-1, 1, len(y))
    U_profile = U_max * np.exp(-x_profile**2 / (2 * width**2))

    # Create 2D velocity fields
    U_jet = np.tile(U_profile[:, np.newaxis], (1, len(x)))
    V_jet = np.zeros_like(U_jet)

    # Build xarray Dataset
    ds_jet = xr.Dataset(
        {
            "u": (["y", "x"], U_jet, {"long_name": "zonal current velocity", "units": "m/s"}),
            "v": (["y", "x"], V_jet, {"long_name": "meridional current velocity", "units": "m/s"}),
        },
        coords={"x": x, "y": y}
    )
    
    return ds_jet