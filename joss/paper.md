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
  - name:Graduate School of Oceanography, University of Rhode Island, Narragansett, RI, USA
    index: 3
date: 7 May 2025
bibliography: paper.bib
---
# Summary
Ocean surface gravity waves are an important component of air-sea interaction, influencing energy, momentum, and gas exchanges across the ocean-atmosphere interface. In specific applications such as refraction by ocean currents or bathymetry, ray tracing provides a computationally efficient way to gain insight into wave propagation. In this paper, we introduce `Mantaray`, an open-source software package implemented in Rust, with a Python interface, that solves the ray equations for ocean surface gravity waves. Mantaray is designed for performance, robustness, and ease of use. The package is modular to facilitate further development and can currently be applied to both idealized and realistic wave propagation problems.

# Statement of need
Ray tracing is a long-standing method for approximating wave propagation across a wide range of disciplines, including optics, seismology, and oceanography, providing a simple framework for studying the evolution of waves  in spatially varying media.

For ocean surface gravity waves, ray-based approaches have been used to study refraction by mesoscale currents, changes in bathymetry, and statistical effects such as directional diffusion of wave action. 

Ray tracing has been widely used in surface wave studies, but the software implementations are often not shared or are written in low-level languages such as Fortran or C, which can be difficult to maintain and integrate into modern workflows. More recently, open-source Python tools—such as the one by Halsne et al. (2023)—have improved accessibility and reproducibility. Mantaray complements these efforts by providing a ray tracing solution built in Rust, a modern  programming language that combines memory safety with execution speed. This choice balances the ease-of-use associated with Python and the computational efficiency of Fortran or C, filling a gap for users who need robust, high-performance ray tracing within a user-friendly environment.

While Rust is still relatively new in the scientific software ecosystem, especially in oceanography, the development of Mantaray illustrates its potential for broader adoption in geoscientific computing. Our package aims to help establish Rust as a top-of-mind language for developing efficient, modern scientific software.

# Key Features
Mantaray is composed of two primary layers:

1. Core Engine (Rust): Implements the numerical integration of the ray equations considering stationary currents.

$$\dot \mathbf{x} =  \mathbf{c_g} + \mathbf{U}$$
$$\dot \mathbf{k} =  -\frac{\nabla \sigma(k)} - \frac{\mathbf{\nabla} \mathbf{k \cdot U}(x, y}}$$



2. Python Interface: Provides a high-level API for initializing simulations, supplying input fields, and running ray integrations. 

# Acknowledgements
ABVB 

# References
