use std::path::Path;

use netcdf3::FileReader;

use super::BathymetryData;
use crate::{error::Error, interpolator};

/// A struct that stores a netcdf3 dataset named test_bathy_3.nc with
/// methods to access, find nearest values, interpolate, and return depth.
///
/// # Note
/// Currently, the methods do not know the difference between an out of
/// bounds point and a just inside the bounds point. The nearest to each of
/// these will be on the edge, so it will return None for both. In the
/// future, the methods should be able to distinguish these two cases.
///
/// In this struct, None is used when the current function will not panic,
/// but the value is not useful to the other structs. Error is used when the
/// function would panic, so instead, it returns an error.
pub(crate) struct CartesianFile {
    variables: (Vec<f32>, Vec<f32>, Vec<f32>),
}

impl BathymetryData for CartesianFile {
    /// Returns the interpolated depth for the given x, y coordinate.
    ///
    /// # Arguments
    /// `x` : `&f32`
    /// - x coordinate
    ///
    /// `y` : `&f32`
    /// - y coordinate
    ///
    /// # Returns
    /// `Result<f32, Error>`
    /// - `Ok(f32)` : depth at the point
    /// - `Err(Error)` : error during execution of `get_depth`.
    ///
    /// # Errors
    /// - `Error::CornerOutOfBounds` : this error is returned when the
    ///   `four_corners` method errors.
    /// - `Error::IndexOutOfBounds` : this error is returned when the `x` or
    /// `y` input give an out of bounds output during the `interpolate`
    /// method.
    /// - `Error::InvalidArgument` : this error is returned from
    ///   `interpolator::bilinear` due to incorrect argument passed.
    /// - `Error::NoNearestPoint` : The target point was either outside the
    /// domain or closest to the edge of the domain.
    fn get_depth(&self, x: &f32, y: &f32) -> Result<f32, Error> {
        if x.is_nan() || y.is_nan() {
            return Ok(f32::NAN);
        }

        let nearest_pt = match self.nearest_point(x, y) {
            Some(p) => p,
            None => return Err(Error::NoNearestPoint),
        };
        let edge_points = match self.four_corners(&nearest_pt.0, &nearest_pt.1) {
            Some(p) => p,
            None => return Err(Error::CornersOutOfBounds),
        };
        let depth = self.interpolate(&edge_points, &(*x, *y))?;
        Ok(depth)
    }

    /// Returns the interpolated depth and gradient for the given x, y coordinate.
    ///
    /// # Arguments
    /// `x` : `&f32`
    /// - x coordinate
    ///
    /// `y` : `&f32`
    /// - y coordinate
    ///
    /// # Returns
    /// `Result<(f32, (f32, f32)), Error>`
    /// - `Ok((f32, (f32, f32)))` : (h, (dhdx, dhdy)), the depth and gradient at the point
    /// - `Err(Error)` : error during execution of `get_depth`.
    ///
    /// # Errors
    /// - `Error::IndexOutOfBounds` : this error is returned when the
    /// `x` or `y` input give an out of bounds output.
    /// - `Error::InvalidArgument` : this error is returned from
    ///   `interpolator::bilinear` due to incorrect argument passed.
    fn get_depth_and_gradient(&self, x: &f32, y: &f32) -> Result<(f32, (f32, f32)), Error> {
        if x.is_nan() || y.is_nan() {
            return Ok((f32::NAN, (f32::NAN, f32::NAN)));
        }

        // find nearest and surrounding edge points
        let nearest_pt = match self.nearest_point(x, y) {
            Some(p) => p,
            None => return Err(Error::IndexOutOfBounds), // TODO: for none should it error or NAN or something else?
        };
        let edge_points = match self.four_corners(&nearest_pt.0, &nearest_pt.1) {
            Some(p) => p,
            None => return Err(Error::IndexOutOfBounds), // TODO: same as above comment
        };

        // interpolate the depth
        let depth = self.interpolate(&edge_points, &(*x, *y))?;

        // get the gradient

        // NOTE: the gradient assumes that the depth is linear in both the x
        // and y directions, and since bilinear interpolation is used to
        // interpolate the depth at any given point, this is a good
        // approximation.
        let x_space = self.variables.0[1] - self.variables.0[0];
        let y_space = self.variables.1[1] - self.variables.1[0];

        let x_grad = (self.depth_from_arr(&edge_points[3].0, &edge_points[3].1)?
            - self.depth_from_arr(&edge_points[1].0, &edge_points[1].1)?)
            / (2.0 * x_space);

        let y_grad = (self.depth_from_arr(&edge_points[2].0, &edge_points[2].1)?
            - self.depth_from_arr(&edge_points[0].0, &edge_points[0].1)?)
            / (2.0 * y_space);

        Ok((depth, (x_grad, y_grad)))
    }
}

impl CartesianFile {
    #[allow(dead_code)]
    /// Construct CartesianFile
    ///
    /// # Arguments
    /// `path` : `&Path`
    /// - a path to the location of the netcdf3 file
    ///
    /// # Returns
    /// `Self` : an initialized CartesianFile
    ///
    /// # Panics
    /// `new` will panic if the data type is invalid or if any of the names are
    /// invalid. But this should never panic for test_bathy_3.nc, unless the
    /// path to this file is incorrect.
    ///
    /// # Note
    /// in the future, be able to check attributes and verify that the file is
    /// correct.
    pub(crate) fn new(path: &Path) -> Self {
        let mut data = FileReader::open(path).unwrap();

        let variables = (
            data.read_var_f32("x").unwrap(),
            data.read_var_f32("y").unwrap(),
            data.read_var_f64("depth")
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
        );
        CartesianFile { variables }
    }

    /// Find nearest point
    ///
    /// # Arguments
    /// `target` : `f32`
    /// - the value to find
    ///
    /// `arr` : `&Vec<f32>`
    /// - the array that will be used when searching for the closest value.
    ///
    /// # Returns
    /// `usize`: index of closest value
    fn nearest(&self, target: &f32, arr: &[f32]) -> usize {
        let mut closest_index = 0;
        let mut closest_distance = (target - arr[0]).abs();

        for (i, value) in arr.iter().enumerate().skip(1) {
            let distance = (target - value).abs();

            if distance < closest_distance {
                closest_index = i;
                closest_distance = distance;
            }
        }
        closest_index
    }

    /// Returns the nearest x index, y index point to given x, y coordinate
    ///
    /// # Arguments
    /// `x`: `&f32`
    /// - x coordinate
    ///
    /// `y`: `&f32`
    /// - y coordinate
    ///
    /// # Returns
    /// `Option<(usize, usize)>`
    /// - `Some((usize, usize))` : the nearest point to the given `x`, `y`
    ///   input.
    /// - `None` : based on the calculated nearest point, the given `x` and
    ///   `y` are assumed to be out of bounds, so a value does not exist.
    ///
    /// # Note
    /// This function will never panic, but if given an out of bounds point,
    /// it will return the closest edge. To attempt to fix this problem,
    /// `nearest_point` should return `None` on points that are on the edges.
    fn nearest_point(&self, x: &f32, y: &f32) -> Option<(usize, usize)> {
        let indx = self.nearest(x, &self.variables.0);
        let indy = self.nearest(y, &self.variables.1);

        if indx == 0
            || indy == 0
            || indx >= self.variables.0.len() - 1
            || indy >= self.variables.1.len()
        {
            return None;
        }

        Some((indx, indy))
    }

    /// Get four adjecent points
    ///
    /// # Arguments
    /// `indx` : `&usize`
    /// - index of the x location
    ///
    /// `indy` : `&usize`
    /// - index of the y location
    ///
    /// # Returns
    /// `Option<Vec<(usize, usize)>>`
    /// - `Some(Vec<(usize, usize)>)` : corner points in surrounding given
    ///   `indx` and `indy` in clockwise order.
    /// - `None` : `indx` or `indy` is out of range and no value exists.
    fn four_corners(&self, indx: &usize, indy: &usize) -> Option<Vec<(usize, usize)>> {
        if *indx == 0
            || *indy == 0
            || *indx >= self.variables.0.len() - 1
            || *indy >= self.variables.1.len() - 1
        {
            return None;
        }
        // clockwise in order
        Some(vec![
            (*indx, indy + 1),
            (indx + 1, *indy),
            (*indx, indy - 1),
            (indx - 1, *indy),
        ])
    }

    /// Interpolate the depth using crate::interpolator::bilinear
    ///
    /// First, the index points are converted to the x and y values at those
    /// indexes, then the depth at that index is taken. Finally, these are
    /// used as arguments to `interpolator::bilinear`.
    ///
    /// # Arguments
    /// `points`: `&Vec<(usize, usize)>`
    /// - a vector of defined points in the depth grid
    ///
    /// `target`: `&(f32, f32)`
    /// - interpolate the depth at this point
    ///
    /// # Returns
    /// `Result<f32, Error>`
    /// - `Ok(f32)` : the depth at the target point
    /// - `Err(Error)` : cannot read depths from at coordinates in the
    ///   `points` vector.
    ///
    /// # Errors
    /// - `Error::IndexOutOfBounds` : one or more of the points passed to
    /// `points` is out of bounds.
    /// - `Error::InvalidArgument` : error during execution of
    /// `interpolator::bilinear` due to invalid arguments.
    fn interpolate(&self, points: &[(usize, usize)], target: &(f32, f32)) -> Result<f32, Error> {
        let pts = vec![
            (
                self.variables.0[points[0].0],
                self.variables.1[points[0].1],
                self.depth_from_arr(&points[0].0, &points[0].1)?,
            ),
            (
                self.variables.0[points[1].0],
                self.variables.1[points[1].1],
                self.depth_from_arr(&points[1].0, &points[1].1)?,
            ),
            (
                self.variables.0[points[2].0],
                self.variables.1[points[2].1],
                self.depth_from_arr(&points[2].0, &points[2].1)?,
            ),
            (
                self.variables.0[points[3].0],
                self.variables.1[points[3].1],
                self.depth_from_arr(&points[3].0, &points[3].1)?,
            ),
        ];
        interpolator::bilinear(&pts, target)
    }

    /// Access values in flattened array as you would a 2d array
    ///
    /// # Arguments
    /// `indx` : `usize`
    /// - index of location in x array
    ///
    /// `indy` : `usize`
    /// - index of location in y array
    ///
    /// # Returns
    /// `Result<f32, Error>`
    /// - `Ok(f32)` : depth
    /// - `Err(Error::IndexOutOfBounds)` : the combined index (x_length *
    ///   indy + indx) is out of bounds of the depth array.
    ///
    /// # Errors
    /// `Err(Error::IndexOutOfBounds)` : this error is returned when `indx`
    /// and `indy` produce a value outside of the depth array.
    fn depth_from_arr(&self, indx: &usize, indy: &usize) -> Result<f32, Error> {
        let index = self.variables.0.len() * indy + indx;
        if index >= self.variables.2.len() {
            return Err(Error::IndexOutOfBounds);
        }
        Ok(self.variables.2[index])
    }
}

#[cfg(test)]
mod test_cartesian_file {

    use crate::{
        bathymetry::{cartesian::CartesianFile, BathymetryData},
        error::Error,
    };
    use std::path::Path;

    /// Create a test file like the test_bathy.nc
    fn create_file(path: &Path, x_len: usize, y_len: usize, x_step: f32, y_step: f32) {
        let x_data: Vec<f32> = (0..x_len).map(|x| x as f32 * x_step).collect();
        let y_data: Vec<f32> = (0..y_len).map(|y| y as f32 * y_step).collect();

        let mut depth_data: Vec<f64> = Vec::new();
        for y in &y_data {
            for x in &x_data {
                if *x < 25000.0 {
                    if *y < 12500.0 {
                        depth_data.push(20.0);
                    } else {
                        depth_data.push(10.0);
                    }
                } else {
                    if *y < 12500.0 {
                        depth_data.push(5.0);
                    } else {
                        depth_data.push(15.0);
                    }
                }
            }
        }

        // most below copied from the docs
        use netcdf3::{DataSet, FileWriter, Version};
        let y_dim_name: &str = "y";
        let y_var_name: &str = y_dim_name;
        let y_var_len: usize = y_len;

        let x_dim_name: &str = "x";
        let x_var_name: &str = x_dim_name;
        let x_var_len: usize = x_len;

        let depth_var_name: &str = "depth";
        let depth_var_len: usize = depth_data.len();

        // Create the NetCDF-3 definition
        // ------------------------------
        let data_set: DataSet = {
            let mut data_set: DataSet = DataSet::new();
            // Define the dimensions
            data_set.add_fixed_dim(y_dim_name, y_var_len).unwrap();
            data_set.add_fixed_dim(x_dim_name, x_var_len).unwrap();
            // Define the variable
            data_set.add_var_f32(y_var_name, &[y_dim_name]).unwrap();
            data_set.add_var_f32(x_var_name, &[x_var_name]).unwrap();
            data_set
                .add_var_f64(depth_var_name, &[y_dim_name, x_var_name])
                .unwrap();

            data_set
        };

        // ...

        // Create and write the NetCDF-3 file
        // ----------------------------------
        let mut file_writer: FileWriter = FileWriter::open(path).unwrap();
        // Set the NetCDF-3 definition
        file_writer.set_def(&data_set, Version::Classic, 0).unwrap();
        assert_eq!(depth_var_len, x_var_len * y_var_len);
        file_writer.write_var_f32(y_var_name, &y_data[..]).unwrap();
        file_writer.write_var_f32(x_var_name, &x_data[..]).unwrap();
        file_writer
            .write_var_f64(depth_var_name, &depth_data[..])
            .unwrap();
        // file_writer.close().unwrap();
        // end of copied from docs
    }

    /// Create a test file like the test_bathy.nc
    fn create_fn_file(
        path: &Path,
        x_len: usize,
        y_len: usize,
        x_step: f32,
        y_step: f32,
        depth_fn: fn(usize, usize) -> f64,
    ) {
        let x_data: Vec<f32> = (0..x_len).map(|x| x as f32 * x_step).collect();
        let y_data: Vec<f32> = (0..y_len).map(|y| y as f32 * y_step).collect();

        let depth_data: Vec<f64> = (0..x_len * y_len)
            .map(|p| depth_fn(p % x_len, p / y_len))
            .collect();

        // most below copied from the docs
        use netcdf3::{DataSet, FileWriter, Version};
        let y_dim_name: &str = "y";
        let y_var_name: &str = y_dim_name;
        let y_var_len: usize = y_len;

        let x_dim_name: &str = "x";
        let x_var_name: &str = x_dim_name;
        let x_var_len: usize = x_len;

        let depth_var_name: &str = "depth";
        let depth_var_len: usize = depth_data.len();

        // Create the NetCDF-3 definition
        // ------------------------------
        let data_set: DataSet = {
            let mut data_set: DataSet = DataSet::new();
            // Define the dimensions
            data_set.add_fixed_dim(y_dim_name, y_var_len).unwrap();
            data_set.add_fixed_dim(x_dim_name, x_var_len).unwrap();
            // Define the variable
            data_set.add_var_f32(y_var_name, &[y_dim_name]).unwrap();
            data_set.add_var_f32(x_var_name, &[x_var_name]).unwrap();
            data_set
                .add_var_f64(depth_var_name, &[y_dim_name, x_var_name])
                .unwrap();

            data_set
        };

        // ...

        // Create and write the NetCDF-3 file
        // ----------------------------------
        let mut file_writer: FileWriter = FileWriter::open(path).unwrap();
        // Set the NetCDF-3 definition
        file_writer.set_def(&data_set, Version::Classic, 0).unwrap();
        assert_eq!(depth_var_len, x_var_len * y_var_len);
        file_writer.write_var_f32(y_var_name, &y_data[..]).unwrap();
        file_writer.write_var_f32(x_var_name, &x_data[..]).unwrap();
        file_writer
            .write_var_f64(depth_var_name, &depth_data[..])
            .unwrap();
        // file_writer.close().unwrap();
        // end of copied from docs
    }

    #[test]
    /// This test checks that the file was created.
    fn test_create_test_bathy() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy1.nc")).unwrap();

        create_file(lockfile.path(), 101, 51, 500.0, 500.0);

        assert_eq!(lockfile.path(), Path::new("tmp_bathy1.nc"))
    }

    #[test]
    // test accessing and viewing variables
    fn test_vars() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy2.nc")).unwrap();

        create_file(lockfile.path(), 101, 51, 500.0, 500.0);

        let data = CartesianFile::new(Path::new(lockfile.path()));
        assert!((data.variables.0[10] - 5000.0).abs() < f32::EPSILON)
    }

    #[test]
    // test the and view the nearest function
    fn test_get_nearest() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy3.nc")).unwrap();

        create_file(lockfile.path(), 101, 51, 500.0, 500.0);

        let data = CartesianFile::new(Path::new(lockfile.path()));
        assert!(data.nearest(&5499.0, &data.variables.0) == 11);
    }

    #[test]
    // check the output from four_corners function
    fn test_get_corners() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy4.nc")).unwrap();

        create_file(lockfile.path(), 101, 51, 500.0, 500.0);

        let data = CartesianFile::new(Path::new(lockfile.path()));
        let corners = data.four_corners(&10, &10).unwrap();
        assert!(corners[0].0 == 10 && corners[0].1 == 11);
        assert!(corners[1].0 == 11 && corners[1].1 == 10);
        assert!(corners[2].0 == 10 && corners[2].1 == 9);
        assert!(corners[3].0 == 9 && corners[3].1 == 10);
    }

    #[test]
    // check values inside the four quadrants but not on grid points
    fn test_depth() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy5.nc")).unwrap();

        create_file(lockfile.path(), 101, 51, 500.0, 500.0);

        let data = CartesianFile::new(Path::new(lockfile.path()));

        // check to see if depth is the same as above
        let check_depth = vec![
            (10099.0, 5099.0, 20.0),
            (30099.0, 5099.0, 5.0),
            (10099.0, 15099.0, 10.0),
            (30099.0, 15099.0, 15.0),
        ];

        for (x, y, h) in &check_depth {
            let depth = data.get_depth_and_gradient(x, y).unwrap().0;
            assert!(
                (depth - h).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                h,
                depth
            );
        }
    }

    #[test]
    /// tests if an IndexOutOfBounds error is returned when accessing depth that
    /// is out of bounds in the x direction
    fn test_x_out_of_bounds() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy6.nc")).unwrap();

        create_file(lockfile.path(), 101, 51, 500.0, 500.0);

        let data = CartesianFile::new(Path::new(lockfile.path()));
        if let Error::NoNearestPoint = data.get_depth(&-500.1, &500.1).unwrap_err() {
            assert!(true);
        } else {
            assert!(false);
        }
    }

    #[test]
    /// tests if an IndexOutOfBounds error is returned when accessing depth that
    /// is out of bounds in the y direction
    fn test_y_out_of_bounds() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy7.nc")).unwrap();

        create_file(lockfile.path(), 101, 51, 500.0, 500.0);

        let data = CartesianFile::new(Path::new(lockfile.path()));
        if let Error::NoNearestPoint = data.get_depth(&500.1, &-500.1).unwrap_err() {
            assert!(true);
        } else {
            assert!(false);
        }
    }

    #[test]
    // test edge cases and center with different depth points. These are
    // using grid points so that it is easy to verify them as the average of
    // the nearest 4 corners.
    fn test_more_depths() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy8.nc")).unwrap();

        create_file(lockfile.path(), 101, 51, 500.0, 500.0);

        let data = CartesianFile::new(Path::new(lockfile.path()));

        // check to see if depth is the same as above
        let check_depth = vec![
            (23000.0, 20000.0, 10.0),
            (10000.0, 12500.0, 12.5),
            (25000.0, 5000.0, 8.75),
            (40000.0, 12500.0, 12.5),
        ];

        for (x, y, h) in &check_depth {
            let depth = data.get_depth_and_gradient(x, y).unwrap().0;
            assert!(
                (depth - h).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                h,
                depth
            );
        }
    }

    #[test]
    fn test_nan() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy9.nc")).unwrap();

        create_file(lockfile.path(), 101, 51, 500.0, 500.0);

        let data = CartesianFile::new(Path::new(lockfile.path()));

        let nan = f32::NAN;

        assert!(data.get_depth(&nan, &nan).unwrap().is_nan());
        assert!(data.get_depth(&10000.0, &nan).unwrap().is_nan());
        assert!(data.get_depth(&nan, &10000.0).unwrap().is_nan());
    }

    #[test]
    fn test_depth_and_gradient_x() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy10.nc")).unwrap();

        fn depth_fn(indx: usize, _indy: usize) -> f64 {
            indx as f64 * 0.05
        }

        create_fn_file(lockfile.path(), 100, 100, 1.0, 1.0, depth_fn);

        let data = CartesianFile::new(Path::new(lockfile.path()));

        // check to see if depth is the same as above
        let check_depth = vec![(10.0, 30.0, 0.5), (30.0, 10.0, 1.5)];

        for (x, y, h) in &check_depth {
            let depth = data.get_depth_and_gradient(x, y).unwrap().0;
            assert!(
                (depth - h).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                h,
                depth
            );
        }

        // check to see if gradient is the same
        let check_gradient = vec![
            (50.0, 50.0, -0.05, 0.0),
            (14.0, 12.0, -0.05, 0.0),
            (10.0, 80.0, -0.05, 0.0),
        ];

        for (x, y, dhdx, dhdy) in &check_gradient {
            let x_grad = data.get_depth_and_gradient(x, y).unwrap().1 .0;
            let y_grad = data.get_depth_and_gradient(x, y).unwrap().1 .1;
            assert!(
                (x_grad - dhdx).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                dhdx,
                x_grad
            );
            assert!(
                (y_grad - dhdy).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                dhdy,
                y_grad
            );
        }
    }

    #[test]
    fn test_depth_and_gradient_y() {
        // create temporary file
        use lockfile::Lockfile;
        let lockfile = Lockfile::create(Path::new("tmp_bathy11.nc")).unwrap();

        fn depth_fn(_indx: usize, indy: usize) -> f64 {
            indy as f64 * 0.05
        }

        create_fn_file(lockfile.path(), 100, 100, 1.0, 1.0, depth_fn);

        let data = CartesianFile::new(Path::new(lockfile.path()));

        // check to see if depth is the same as above
        let check_depth = vec![(10.0, 30.0, 1.5), (30.0, 10.0, 0.5)];

        for (x, y, h) in &check_depth {
            let depth = data.get_depth_and_gradient(x, y).unwrap().0;
            assert!(
                (depth - h).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                h,
                depth
            );
        }

        // check to see if gradient is the same
        let check_gradient = vec![
            (50.0, 50.0, 0.0, -0.05),
            (14.0, 12.0, 0.0, -0.05),
            (10.0, 80.0, 0.0, -0.05),
        ];

        for (x, y, dhdx, dhdy) in &check_gradient {
            let x_grad = data.get_depth_and_gradient(x, y).unwrap().1 .0;
            let y_grad = data.get_depth_and_gradient(x, y).unwrap().1 .1;
            assert!(
                (x_grad - dhdx).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                dhdx,
                x_grad
            );
            assert!(
                (y_grad - dhdy).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                dhdy,
                y_grad
            );
        }
    }
}
