import xarray as xr
import numpy as np
import numpy.matlib

import time
import math

import matplotlib.pyplot as plt


creator_name = 'g.marechal'

g = 9.8
def frequency_from_wavelength(lambda_p, depth):
    """
    Purpose: 
    ---------
    f from full dispersion relationship
    """
    
    k = 2*np.pi/lambda_p # the wavenumber
    f = (g*k * np.tanh(g*depth))/(2*np.pi)
    period = 1/f
    return f

    
def group_velocity(k, f, depth):
    """
    Purpose: 
    ---------
    Compute Cg from the full dispersion relationship
    """
    sigma = 2*np.pi*f
    phase_velocity = sigma/k
    cg = phase_velocity * (1/2 + (k*depth)/(np.sinh(2*k*depth)))
    return cg


def find_closest_pixel(x_line_param, y_line_param, X, Y):
    """
    Purpose: 
    ---------
    Function to find the closest grid point to a line point (x_line_param, y_line_param) 
    """


    distances = np.sqrt((X - x_line_param)**2 + (Y - y_line_param)**2)
    idx = np.unravel_index(np.argmin(distances), distances.shape)  # Find the index of the closest point

    return X[idx], Y[idx]


def find_closest_pixel_depth(x_line_param, y_line_param, X, Y, depth):
    distances = np.sqrt((X - x_line_param)**2 + (Y - y_line_param)**2)
    idx = np.unravel_index(np.argmin(distances), distances.shape)  # Find the index of the closest point

    return X[idx], Y[idx], depth[idx]

    
def start_line_ray_tracing(x_line, y_line, x_grid, y_grid, num_points):
    
    """
    Purpose: 
    ---------
    The starting line for ray tracing computation
    """
    
    x_line_param = np.linspace(x_line[0], x_line[1], num_points)
    y_line_param = np.linspace(y_line[0], y_line[1], num_points)

    # Loop through the points on the line and find the closest pixels
    closest_x = []
    closest_y = []
    closest_depth = []
    XX, YY = np.meshgrid(x_grid, y_grid)
    
    for x, y in zip(x_line_param, y_line_param):
        x_closest, y_closest = find_closest_pixel(x, y, XX, YY)
        closest_x.append(float(x_closest))
        closest_y.append(float(y_closest))
        
    return closest_x, closest_y

def point_to_line_distance(x, y, x1, y1, x2, y2):
    # Line equation: (y2 - y1)(x - x1) = (x2 - x1)(y - y1)
    num = abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1)
    denom = np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
    return num / denom



    
def create_current_forcing(x, y, ucur, vcur, output_file_cur, output_file_depth, depth =  4e3):
    """
    Purpose: 
    ---------
    create Xarray DataArray forcing for Ray tracing.
    Crate a Netcdf file with the zonal and meridional currents and constant depth
    """
    # --- Current
    ds = xr.Dataset(
    data_vars={
        "x": (["x"], x),
        "y": (["y"], y),
        "u": (["y", "x"], ucur, {"long_name": "u", "units": "m/s"}),
        "v": (["y", "x"], vcur, {"long_name": "v", "units": "m/s"}),
    })

    ds.to_netcdf(output_file_cur, format="NETCDF3_CLASSIC")

    # --- Depth
    depth = depth*np.ones((len(y), len(x)))
    
    ds = xr.Dataset(
        data_vars={
            "x": (["x"], x),
            "y": (["y"], y),
            "depth": (["y", "x"], depth, {"long_name": "depth", "units": "m"}),
        })
    
    ds.to_netcdf(output_file_depth, format="NETCDF3_CLASSIC")
    
def create_depth_forcing(x, y, depth, output_file_cur, output_file_depth, ucur = 0, vcur = 0):
    """
    Purpose:
    ---------
    create Xarray DataArray forcing for Ray tracing. Crate a Netcdf file with the depth and constant currents (= 0)
    """

        # --- Depth
    ds = xr.Dataset(
        data_vars={
            "x": (["x"], x),
            "y": (["y"], y),
            "depth": (["y", "x"], depth, {"long_name": "depth", "units": "m"}),
        })
    
    ds.to_netcdf(output_file_depth, format="NETCDF3_CLASSIC")
    
    # --- Current
    ucur = ucur*np.zeros((len(y), len(x)))
    vcur = vcur*np.zeros((len(y), len(x)))

    ds = xr.Dataset(
    data_vars={
        "x": (["x"], x),
        "y": (["y"], y),
        "u": (["y", "x"], ucur, {"long_name": "u", "units": "m/s"}),
        "v": (["y", "x"], vcur, {"long_name": "v", "units": "m/s"}),
    })

    ds.to_netcdf(output_file_cur, format="NETCDF3_CLASSIC")

def decimal_coords_to_meters(delta_lon, delta_lat, mean_latitude):
    """
    Purpose:
    ---------
    Convert decimal longitude and latitude to meters.
    
    Inputs:
    ---------
    - delta_lon: Change in longitude in decimal degrees.
    - delta_lat: Change in latitude in decimal degrees.
    - mean_latitude: Mean latitude in decimal degrees.
   Outputs:
   ---------
    - (meters_x, meters_y): Distances in meters (x: east-west, y: north-south).
    """
    # Conversion constants
    meters_per_degree_lat = 111320  # Approximate meters per degree latitude
    meters_per_degree_lon = 111320 * np.cos(np.radians(mean_latitude))  # Adjust for latitude
    
    # Convert decimal degrees to meters
    meters_y = delta_lat * meters_per_degree_lat
    meters_x = delta_lon * meters_per_degree_lon
    
    return meters_x, meters_y


