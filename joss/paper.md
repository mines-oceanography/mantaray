TODO: update the joss required yaml section

---
title: TODO
tags: TODO
authors:
    - name: Bryce Irving
      orcid: 0009-0004-2309-9522
      affiliation: 1
    - name: Guilherme P. Castelao
      orcid: 0000-0002-6765-0708
      affiliation: 2
    - name: Bia Villas Boas
      orcid: TODO
      affiliation: 1
affiliations:
    - name: Colorado School of Mines
      index: 1
    - name: TODO
      index: 2
date: 28 May 2024 (TODO this might need updating?)
bibliography: paper.bib
---

# Summary

Ray-tracing is a powerful technique used in computer graphics and scientific simulations to model the propagation of waves. When it comes to ocean surface gravity waves, ray-tracing can provide valuable insights into wave propagation, including their interaction with ocean currents and the ocean bathymetry. Such insights are crucial for better understanding ocean wave physics and their role in climate, as waves are a major player in heat transfer and the exchange of gases between the ocean and the atmosphere.
In the present work, we introduce a novel package to calculate the path of a wave propagating through the ocean. We developed the code in Rust for memory safety, high performance, robustness, and simple testing. The ray tracing uses the Runge-Kutta 4th order and bilinear interpolation methods to reduce integration error. The package includes supporting Python files to visualize the results. The components are tested individually, and the overall output is checked against known idealized cases. With these methods, our Rust crate can trace the propagation of a wave through a variable depth represented in cartesian coordinates and plot the results. These results are significant because accurately tracing the propagation of ocean waves with an efficient language will increase performance making it seamless to run large simulation ensembles. There are many reasonable opportunities for improvement in the future. The accuracy and realism of the program will improve by tracing bundles of multiple rays and accounting for the interactions with ocean currents. Additionally, the computational performance will improve by parallelizing the ray tracing.

# Statement of need

# References
