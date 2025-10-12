<!-- start elevator-pitch -->

[![Rust checks](https://github.com/mines-oceanography/ray_tracing/actions/workflows/ci.yml/badge.svg)](https://github.com/mines-oceanography/ray_tracing/actions/workflows/ci.yml)

# Ray Tracing

A library for surface gravity waves ray tracing.

![Demo](https://github.com/mines-oceanography/mantaray/blob/main/notebooks/canonical_examples/demo_animation.gif)

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

2. Clone the repo
```
git clone git@github.com:mines-oceanography/mantaray.git
cd mantaray
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

*   **Reporting bugs:**  If you find a bug, please open an [issue](https://github.com/mines-oceanography/mantaray/issues) with the `bug` label and provide as much detail as possible, including steps to reproduce the issue.
*   **Suggesting features:**  Have an idea for a new feature or improvement? Open an [issue](https://github.com/mines-oceanography/mantaray/issues) with the `enhancement` label and describe your suggestion.
*   **Submitting code changes:**  We welcome code contributions!  Please follow the Pull Request Guidelines below.
*   **Improving documentation:**  Clear and concise documentation is essential. If you find areas where the documentation can be improved, please submit an [issue](https://github.com/mines-oceanography/mantaray/issues) with the `documentation` label.

### Pull Request Guidelines

Before submitting a pull request, please make sure it meets these guidelines:

1.  **Tests:**  All pull requests should include unit tests that cover the changes.
2.  **Documentation:**  If your pull request adds or modifies functionality, please update the documentation accordingly.
3.  **CI:**  Your pull request must pass all existing continuous integration checks.
4.  **Single Functionality:**  Each pull request should ideally address a single, well-defined functionality.  If your changes are more extensive, please consider breaking them down into multiple, smaller pull requests.

### Getting Help

If you have questions or need help getting started, please open an issue with the `question` label.  We'll do our best to assist you.

<!-- end elevator-pitch -->
