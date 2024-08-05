import numpy as np
import xarray as xr

slope = 0.1
L = 10000
dx = 10
r0 = 500 # Island diameter = 1000m
r_out = 8e3
max_depth = r_out*slope
h0 = - r0*slope

x = range(-L,L+1,dx)
x = np.array(x).astype(np.float32)
y = x
Y, X = np.meshgrid(x,y)
R = (X**2 + Y**2)**0.5
h = max_depth*np.ones(R.shape) + h0
ind = R <= r_out
h[ind] = R[ind] * slope + h0
h = np.array(h).astype(np.float64)

ds = xr.Dataset(
    data_vars={
        "x": (["x"], x),
        "y": (["y"], y),
        "depth": (["x", "y"], h, {"long_name": "Depth", "units": "m"}),
    })

ds.to_netcdf("island.nc", format="NETCDF3_CLASSIC")
