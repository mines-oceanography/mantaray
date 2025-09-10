Notebook Gallery
================

These notebooks contain examples of various use cases for mantaray.

Each example is in its own folder. Many folders have more than one notebook, and some will require you to run a notebook that generates data first before running the ray tracing notebook. 

To run the example notebooks, follow installation instructions in the readme to install pixi and clone the repo. Then, install the examples environment using `pixi run -e examples develop`. 

After that there are multiple options:

1. Follow our development instructions for using jupyter lab.

2. Start a shell with the environment using the command `pixi shell -e examples`.

3. Find the installation inside the `.pixi` folder and run the examples your own way.

Any additional instructions, such as order of running notebooks or needing data files, will be located in the readme for that example in the `notebooks` folder of the repository.

.. base-gallery::
    :caption: Gallery caption
    :tooltip:

    idealized_fields
    canonical_ray_tracing
    homepage_animation
    data_generation
    snells_law_verification
    agulhas_example
    nazare_example
    tutorial