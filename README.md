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

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](https://github.com/mines-oceanography/ray_tracing/blob/main/LICENSE-APACHE "Apache License 2.0") or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](https://github.com/mines-oceanography/ray_tracing/blob/main/LICENSE-MIT "MIT License") or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

<!-- end elevator-pitch -->
