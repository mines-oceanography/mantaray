---
title: 'Mantaray: A Rust Package for Ray Tracing Ocean Surface Gravity Waves'
tags:
  - Rust
  - Ocean
  - Waves
authors:
  - name: Bryce Irving
    orcid: 0009-0004-2309-9522
    affiliation: 1
  - name: Guilherme P. Castelao
    orcid: 0000-0002-6765-0708
    affiliation: 2
  - name: Colin Beyers
    orcid: 0009-0004-8312-6158
    affiliation: 1
  - name: James Clemson
    orcid: 0009-0000-4329-6575
    affiliation: 1
  - name: Jackson Krieger
    orcid: 0009-0006-3693-8887
    affiliation: 1
  - name: Gwendal Marechal
    orcid: 0000-0003-0378-5694
    affiliation: 1
  - name: Nicholas Pizzo
    orcid: 0000-0001-9570-4200
    affiliation: 3
  - name: Bia Villas Bôas
    orcid: 0000-0001-6767-6556
    affiliation: 1
affiliations:
  - name: Colorado School of Mines, Golden, CO, USA
    index: 1
  - name: National Renewable Energy Laboratory, Golden, CO, USA
    index: 2
  - name: Graduate School of Oceanography, University of Rhode Island, Narragansett, RI, USA
    index: 3
date: 7 May 2025
bibliography: paper.bib
---
# Summary
Ocean surface gravity waves are an important component of air-sea interaction, influencing energy, momentum, and gas exchanges across the ocean-atmosphere interface. In specific applications such as refraction by ocean currents or bathymetry, ray tracing provides a computationally efficient way to gain insight into wave propagation. In this paper, we introduce `Mantaray`, an open-source software package implemented in Rust, with a Python interface, that solves the ray equations for ocean surface gravity waves. Mantaray is designed for performance, robustness, and ease of use. The package is modular to facilitate further development and can currently be applied to both idealized and realistic wave propagation problems \ref{fig:examples}.

![Examples of ray tracing performed using `Mantaray` for waves in deep water interacting with a zonal jet (top left) and a mesoscale eddy (bottom left), and for waves encountaring varyring bathymethy on a linear beach (tipo right) and Gaussian island (bottom right) \label{fig:examples}](idealized_showcase.png){ width=100% }

# Statement of need
Ray tracing is a long-standing method for approximating wave propagation across a wide range of disciplines, including optics, seismology, and oceanography, providing a simple framework for studying the evolution of waves  in spatially varying media. For ocean surface gravity waves, ray-based approaches have been used to study refraction by mesoscale currents (e.g., @mapp1985wave, @romero2017observations, @marechal2022variability), changes in bathymetry (e.g., @munk1947refraction, @kukulka2017surface), and statistical effects such as directional diffusion of wave action (e.g., @smit2019swell, @VBY2023). 

Ray tracing has been widely used in surface wave studies, but the software implementations are often not shared or are written in low-level languages such as Fortran or C (e.g., @oreilly2016california), which can be difficult to maintain and integrate into modern workflows. More recently, open-source Python tools—such as the one by @halsne2023ocean—have improved accessibility and reproducibility. Mantaray complements these efforts by providing a ray tracing solution built in Rust, a modern  programming language that combines memory safety and execution speed with out-of-the-box Python integration. This choice balances the ease-of-use associated with Python and the computational efficiency of Fortran or C/C++, filling a gap for users who need robust, high-performance ray tracing within a user-friendly environment.

While Rust is still relatively new in the scientific software ecosystem, especially in oceanography, the development of Mantaray illustrates its potential for broader adoption in geoscientific computing. Our package aims to help establish Rust as a top-of-mind language for developing efficient, modern scientific software.

# Key Features
Mantaray is composed of two primary layers:

1. Core Engine (Rust): Implements the numerical integration of the ray equations considering stationary (no time dependence) currents ${\mathbf U}(x, y)$ in a Cartesian domain.  

	For linear surface gravity waves, with dispersion relationship given by:
	
	$$\sigma = [gk\tanh{(kH(x, y))}]^{1/2},$$
	
	where $\sigma$ is the intrinsic frequency of the waves, $g$ is the gravitational acceleration, $k$ is the wavenumber magnitude, and $H$ is the water depth, the ray equations can be written as:
	
	$${\mathbf {c_g}} = \frac{\partial \sigma}{\partial \mathbf k},$$
	
	$$\dot {\mathbf x} =  {\mathbf {c_g}} + {\mathbf U}(x, y),$$
	
	$$\dot {\mathbf k} =  -{\mathbf \nabla} \sigma -{\mathbf \nabla} \left ( {\mathbf k} \cdot {\mathbf U}\right),$$
	
	where  ${\mathbf c_g}$ is the group velocity,  ${\mathbf k} = (k_x, k_y)$ is the wavenumber vector, and ${\mathbf x} = (x, y)$ is the wave position vector.
	
	`Mantaray` integrates the ray equations using a 4th-order Runge-Kutta scheme from the `ode_solvers` crate, with bilinear interpolation for spatial fields such as bathymetry and surface currents.


2. Python Interface: Provides a high-level API for initializing simulations, supplying input fields, and running ray integrations. The current version of the package supports:

	- Cartesian domains with arbitrary bathymetry and current fields input as NetCDF3 files
	    
	- Configurable integration parameters (step size, duration of the integration)
	    
	- Output of ray paths as Xarray Datasets for easy visualization and diagnostics

`Mantaray` has two main functionalities: `single_ray`, for tracing an individual ray, and `ray_tracing`, for tracing a collection of rays.

*Example:* The following example illustrates the use of the `single_ray` functionality for tracing a wave that is initially propagating from left to right with a wavelength of 100 m. Note that `bathymetry` and `current` are strings with the path to the respective forcing fields.

 ```python
 import numpy as np
 import mantaray

# Define initial conditions 
k0 = 2*np.pi/100 # initial wavenumber magnitude
theta0 = 0 # initial direction
# Calculates wavenumber components
kx0 = k0*np.cos(phi0)
ky0 = k0*np.sin(phi0)

# Define initial position
x0 = 0
y0 = 500

# Define integration parameters
duration = 1000 # duration in seconds
timestep = 0.1 #timestep in seconds

# Performs integration
ray_path = mantaray.single_ray(x0, y0, kx0, ky0, duration, timestep, bathymetry, current)
 ```

Example:* The  `ray_tracing` functionality works similarly, but it takes a collection of initial conditions as `numpy` arrays. In the case below, we are propagating four identical rays, with different initial positions.

 ```python
 import numpy as np
 import mantaray

# Define initial conditions 
k0 = 2*np.pi/100 # initial wavenumber magnitude
theta0 = 0 # initial direction
# Calculates wavenumber components
kx0 = k0*np.cos(phi0)*np.ones(4)
ky0 = k0*np.sin(phi0)*np.ones(4)

# Define initial position
x0 = np.array([0, 0, 0, 0])
y0 = np.array([100, 300, 500, 700])

# Define integration parameters
duration = 1000 # duration in seconds
timestep = 0.1 #timestep in seconds

# Performs integration
ray_path = mantaray.ray_tracing(x0, y0, kx0, ky0, duration, timestep, bathymetry, current)
 ```


# Acknowledgements
We thank the Colorado School of Mines Summer Undergraduate Research Fellowship (SURF) and the Mines Undergraduate Research Fellowship (MURF) for partially supporting undergraduate students BI, JC, and JK. BVB was supported by the ONR MURI award N00014-24-1-2554, and NASA award 80NSSC24K1640 through the SWOT Science Team. CB was supported by NASA awards 80NSSC23K0979 through the International Ocean Vector Winds Science Team and 80GSFC24CA067 through the ODYSEA science team as part of the Earth System Explorer program. JK and GM were supported by NASA award 80NSSC24K0411 through the S-MODE Science Team. All co-authors are thankful to Luiz Irber for helpful recommendations throughout the development of this package.

# References
