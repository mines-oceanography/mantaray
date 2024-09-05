[![Rust checks](https://github.com/mines-oceanography/ray_tracing/actions/workflows/ci.yml/badge.svg)](https://github.com/mines-oceanography/ray_tracing/actions/workflows/ci.yml)

# Ray Tracing

A library for surface gravity waves ray tracing.

## Development
1. Install [Pixi](https://pixi.sh/latest/)

2. Clone the repo
```
git clone git@github.com:mines-oceanography/ray_tracing.git
cd ray_tracing
```

3. Build Python
```
pixi run maturin develop
```
This will take about 20 to 30 minutes (at least on windows 10).

4. Import `mantaray`
```python
from mantaray.core import single_ray, ray_tracing
```
Documentation for these functions are located in [core.py](python/mantaray/core.py).

5. Run Python file using Pixi

```
pixi run python path_to_file.py
```

6. To test Python library run:

```
pixi run -e dev pytest
```



