<!-- start elevator-pitch -->

[![Rust checks](https://github.com/mines-oceanography/ray_tracing/actions/workflows/ci.yml/badge.svg)](https://github.com/mines-oceanography/ray_tracing/actions/workflows/ci.yml)

# Ray Tracing

A library for surface gravity waves ray tracing.

## Development
### Installation
1. Install [Pixi](https://pixi.sh/latest/)

2. Clone the repo
```
git clone git@github.com:mines-oceanography/ray_tracing.git
cd ray_tracing
```

3. Build Python
```
pixi run develop
```
This will take about 20 to 30 minutes (at least for first time compiling on windows 10).

### Usage
At the top of your python file, you will need to include the following import line:
```python
from mantaray.core import single_ray, ray_tracing
```
Documentation for these functions are located in [core.py](#api).

#### Run Python file

```
pixi run python path_to_file.py
```

### Using Jupyter Lab
1. Develop the code for the `jupyterlab` environment
```
pixi run -e jupyterlab develop
```
2. Open Jupyter Lab using the `jupyterlab` environment
```
pixi run -e jupyterlab jupyter lab
```

### To test Python library run:

```
pixi run -e test pytest
```
<!-- end elevator-pitch -->


