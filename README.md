
[![Pixi Badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)
[![Rust checks](https://github.com/mines-oceanography/ray_tracing/actions/workflows/ci.yml/badge.svg)](https://github.com/mines-oceanography/ray_tracing/actions/workflows/ci.yml)

[![Docs](https://img.shields.io/badge/docs-examples-orange)](https://mines-oceanography.github.io/mantaray)
[![Crates.io](https://img.shields.io/crates/v/mantaray)](https://crates.io/crates/mantaray)
[![PyPI version](https://img.shields.io/pypi/v/mantaray)](https://pypi.org/project/mantaray/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)




# Ray Tracing

A library for surface gravity waves ray tracing.

![Demo](https://github.com/mines-oceanography/mantaray/blob/main/notebooks/canonical_examples/demo_animation.gif)

<!-- examples -->

## Examples
The examples are located in the `notebooks` directory, and each scenario is inside its own subfolder.

To run the example notebooks, follow [installation](#installation) instructions to install pixi and clone the repo. Then, install the examples environment using `pixi run -e examples develop`. 

After that there are multiple options:
- Follow our development [instructions](#using-jupyter-lab) for using jupyter lab.
- Start a shell with the environment using the command `pixi shell -e examples`.
- Find the installation inside the `.pixi` folder and run the examples your own way.

If there are additional instructions, such as needing data files, they will be located in the readme of that example's folder.

## Development

### Installation
1. Install [Pixi](https://pixi.sh/latest/)


2. [Fork Mantaray](https://github.com/mines-oceanography/mantaray/fork)'s repository, by clicking in the 'Fork' button in the top-right corner.

3. Clone your forked repository. Check the green button and choose a protocol.
For instance, if you use SSH you'll see something similar to:
```
git clone git@github.com:<your-username>/mantaray.git
cd mantaray
```

4. Build Python
```
pixi run develop
```
This can take a few minutes the very first time.

### Usage
At the top of your python file, you will need to include the following import line:
```python
from mantaray.core import single_ray, ray_tracing
```
Documentation for these functions are located in [core.py](https://mines-oceanography.github.io/mantaray/api.html).

#### Run Python file

```
pixi run python path_to_file.py
```

### Using Jupyter Lab
1. Develop the code for the `examples` environment
```
pixi run -e examples develop
```
2. Open Jupyter Lab using the `examples` environment
```
pixi run -e examples jupyter lab
```

### To test Python library run:

```
pixi run -e test pytest
```

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](https://github.com/mines-oceanography/mantaray/blob/main/LICENSE-APACHE "Apache License 2.0") or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](https://github.com/mines-oceanography/mantaray/blob/main/LICENSE-MIT "MIT License") or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

We welcome contributions to this project!  Whether you're fixing a bug, adding a new feature, or improving the documentation, your help is greatly appreciated. All contributions should be made through GitHub, by forking the repository, creating a new branch, and submitting a pull request.

### Ways to Contribute

There are many ways to contribute to this project, including:

*   **Reporting bugs:**  If you find a bug, please open an [issue](https://github.com/mines-oceanography/mantaray/issues) and provide as much detail as possible, including steps to reproduce the issue.
*   **Suggesting features:**  Have an idea for a new feature or improvement? Open an [issue](https://github.com/mines-oceanography/mantaray/issues) and describe your suggestion.
*   **Submitting code changes:**  We welcome code contributions! If your change is based on an existing issue, please comment on that issue and let us know you are working on it. Otherwise, if it is something new, create an issue and let us know what you are working on. When ready to submit, please follow the Pull Request Guidelines below.
*   **Improving documentation:**  Clear and concise documentation is essential. If you find areas where the documentation can be improved, please submit an [issue](https://github.com/mines-oceanography/mantaray/issues).

> When you create an issue, we may label it (`bug`, `enhancement`, etc). If we are unsure about what you are requesting, we will ask to clarify, and if you believe another label fits it better, let us know.

### Pull Request Guidelines

Before submitting a pull request, please make sure it meets these guidelines:

1.  **Tests:**  All pull requests should include unit tests that cover the changes.
2.  **Documentation:**  If your pull request adds or modifies functionality, please update the documentation accordingly.
3.  **CI:**  Your pull request must pass all existing continuous integration checks.
4.  **Single Functionality:**  Each pull request should ideally address a single, well-defined functionality.  If your changes are more extensive, please consider breaking them down into multiple, smaller pull requests.

### Getting Help

If you have questions or need help getting started, please open an issue and ask.  We'll do our best to assist you.

<!-- end elevator-pitch -->
