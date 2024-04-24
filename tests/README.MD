# Integration tests

The rust test modules are located in this directory and the data files used are
located in the `data/` folder.

This test module currently contains:
 - Linear Sea Mount tests: this test verifies that the rays curve towards the shore when encountering an island and has files:
    - `linear_sea_mount.rs`: contains a test that verifies the rays above the island bend down and the rays below the island bend up.
    - `data/island_slice.nc`: a small slice from the netcdf file produced by `support/linear_sea_mount.py`.

