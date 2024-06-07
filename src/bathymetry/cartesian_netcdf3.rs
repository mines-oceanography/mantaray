//! Struct used to create and access bathymetry data stored in a netcdf3 file.

use std::path::Path;

use netcdf3::{DataType, FileReader};

use super::BathymetryData;
use crate::{error::Error, interpolator};

/// A struct that stores a netcdf3 dataset with methods to access, find nearest
/// values, interpolate, and return depth and gradient.
///
/// # Note
/// Currently, the methods do not know the difference between an out of bounds
/// point and a point within one grid space from the edge. The nearest to each
/// of these will be on the edge, so it will return None for both. In the
/// future, the methods should be able to distinguish these two cases.
///
/// In this struct, None is used when the function will not panic, but
/// the value is not useful to the other structs. Error is used when the
/// function would panic, so instead, it returns an error.
pub struct CartesianNetcdf3 {
    /// a vector containing the x values from the netcdf3 file
    x: Vec<f32>,
    /// a vector containing the y values from the netcdf3 file
    y: Vec<f32>,
    /// a vector containing the depth values from the netcdf3 file. Note this is
    /// a flattened 2d array and is accessed by the function `depth_from_array`.
    depth: Vec<f32>,
}

impl BathymetryData for CartesianNetcdf3 {
    /// Depth at the inputted (x ,y) point.
    ///
    /// # Arguments
    /// `x` : `&f32`
    /// - x location \[m\]
    ///
    /// `y` : `&f32`
    /// - y location \[m\]
    ///
    /// # Returns
    /// `Result<f32, Error>`
    /// - `Ok(f32)` : depth at the point in meters
    /// - `Err(Error)` : error during execution of `get_depth`.
    ///
    /// # Errors
    /// - `Error::IndexOutOfBounds` : this error is returned when the `x` or `y`
    /// input give an out of bounds output during the `interpolate` method.
    /// - `Error::InvalidArgument` : this error is returned from
    ///   `interpolator::bilinear` due to incorrect argument passed.
    fn depth(&self, x: &f32, y: &f32) -> Result<f32, Error> {
        if x.is_nan() || y.is_nan() {
            return Ok(f32::NAN);
        }

        let (nearest_x, nearest_y) = match self.nearest_point(x, y) {
            Some(point) => point,
            None => return Err(Error::IndexOutOfBounds),
        };
        let corner_points = match self.four_corners(&nearest_x, &nearest_y) {
            Some(point) => point,
            None => return Err(Error::IndexOutOfBounds),
        };
        let depth = self.interpolate(&corner_points, &(*x, *y));

        depth
    }

    /// Depth and gradient at the given (x ,y) coordinate.
    ///
    /// # Arguments
    /// `x` : `&f32`
    /// - x location \[m\]
    ///
    /// `y` : `&f32`
    /// - y location \[m\]
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
    fn depth_and_gradient(&self, x: &f32, y: &f32) -> Result<(f32, (f32, f32)), Error> {
        if x.is_nan() || y.is_nan() {
            return Ok((f32::NAN, (f32::NAN, f32::NAN)));
        }

        let (nearest_x, nearest_y) = match self.nearest_point(x, y) {
            Some(point) => point,
            None => return Err(Error::IndexOutOfBounds),
        };
        let corner_points = match self.four_corners(&nearest_x, &nearest_y) {
            Some(point) => point,
            None => return Err(Error::IndexOutOfBounds),
        };

        // interpolate the depth
        let depth = self.interpolate(&corner_points, &(*x, *y))?;

        // get the gradient

        // NOTE: the gradient assumes that the depth is linear in both the x
        // and y directions, and since bilinear interpolation is used to
        // interpolate the depth at any given point, this is a good
        // approximation.
        let x_space = self.x[1] - self.x[0];
        let y_space = self.y[1] - self.y[0];

        let north_point = &corner_points[0];
        let east_point = &corner_points[1];
        let south_point = &corner_points[2];
        let west_point = &corner_points[3];

        let x_gradient = (self.depth_at_indexes(&east_point.0, &east_point.1)?
            - self.depth_at_indexes(&west_point.0, &west_point.1)?)
            / (2.0 * x_space);

        let y_gradient = (self.depth_at_indexes(&north_point.0, &north_point.1)?
            - self.depth_at_indexes(&south_point.0, &south_point.1)?)
            / (2.0 * y_space);

        Ok((depth, (x_gradient, y_gradient)))
    }
}

impl CartesianNetcdf3 {
    #[allow(dead_code)]
    /// Initialize the CartesianNetCDF3 struct with the data from the netcdf3
    /// file
    ///
    /// # Arguments
    /// `path` : `&Path`
    /// - a path to the location of the netcdf3 file
    ///
    /// `xname` : `&str`
    /// - the name of the x variable in the netcdf3 file
    ///
    /// `yname` : `&str`
    /// - the name of the y variable in the netcdf3 file
    ///
    /// `depth_name` : `&str`
    /// - the name of the depth variable in the netcdf3 file
    ///
    /// # Returns
    /// `Self` : an initialized CartesianNetCDF3 struct
    ///
    /// # Panics
    /// `new` will panic if the data type is invalid or if any of the names are
    /// invalid.
    ///
    /// # Note
    /// in the future, be able to check attributes and verify that the file is
    /// correct.
    pub fn open(path: &Path, xname: &str, yname: &str, depth_name: &str) -> Self {
        let mut data = FileReader::open(path).unwrap();

        let x = data.read_var(xname).unwrap();
        let x = match x.data_type() {
            DataType::I16 => x
                .get_i16_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::I8 => x
                .get_i8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::U8 => x
                .get_u8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::I32 => x
                .get_i32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::F32 => x.get_f32_into().unwrap(),
            DataType::F64 => x
                .get_f64_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
        };

        let y = data.read_var(yname).unwrap();
        let y = match y.data_type() {
            DataType::I16 => y
                .get_i16_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::I8 => y
                .get_i8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::U8 => y
                .get_u8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::I32 => y
                .get_i32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::F32 => y.get_f32_into().unwrap(),
            DataType::F64 => y
                .get_f64_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
        };

        let depth = data.read_var(depth_name).unwrap();
        let depth = match depth.data_type() {
            DataType::I16 => depth
                .get_i16_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::I8 => depth
                .get_i8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::U8 => depth
                .get_u8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::I32 => depth
                .get_i32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
            DataType::F32 => depth.get_f32_into().unwrap(),
            DataType::F64 => depth
                .get_f64_into()
                .unwrap()
                .iter()
                .map(|x| *x as f32)
                .collect(),
        };

        CartesianNetcdf3 {
            x,
            y,
            depth,
        }
    }

    /// Find the index of the closest value to the target in the array
    ///
    /// # Arguments
    /// `target` : `&f32`
    /// - the value to find
    ///
    /// `arr` : `&[f32]`
    /// - the array that will be used when searching for the closest value.
    ///
    /// # Returns
    /// `usize`: index of closest value
    ///
    /// # Note
    /// This function uses binary search, but requires the array to be sorted.
    fn nearest(&self, target: &f32, array: &[f32]) -> usize {
        let mut left = 0;
        let mut right = array.len() - 1;
        let mut closest_index = 0;
        let mut closest_distance = f32::MAX;

        // edge cases
        if *target <= array[left] {
            return left;
        }
        if *target >= array[right] {
            return right;
        }

        // binary search
        while left <= right {
            let mid = (left + right) / 2;
            let distance = (target - array[mid]).abs();

            if distance < closest_distance {
                closest_index = mid;
                closest_distance = distance;
            }

            if array[mid] == *target {
                return mid;
            }

            if array[mid] > *target {
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }

        closest_index
    }

    /// Returns the nearest (xindex, yindex) point to given (x ,y) point
    ///
    /// # Arguments
    /// `x`: `&f32`
    /// - x location in meters
    ///
    /// `y`: `&f32`
    /// - y location in meters
    ///
    /// # Returns
    /// `Option<(usize, usize)>`
    /// - `Some((usize, usize))` : the indices of the nearest point to the given
    ///   `x`, `y` input.
    /// - `None` : based on the calculated nearest point, the given `x` and `y`
    ///   are assumed to be out of bounds, so a value does not exist.
    ///
    /// # Note
    /// This function will return `None` on points that are on the edges. This
    /// causes a small bug requiring the user to initialize the ray at least
    /// half a grid space away from the edge.
    fn nearest_point(&self, x: &f32, y: &f32) -> Option<(usize, usize)> {
        let xindex = self.nearest(x, &self.x);
        let yindex = self.nearest(y, &self.y);

        if xindex == 0 || yindex == 0 || xindex >= self.x.len() - 1 || yindex >= self.y.len() {
            return None;
        }

        Some((xindex, yindex))
    }

    /// Get four adjacent points
    ///
    /// # Arguments
    /// `xindex` : `&usize`
    /// - index of the x location
    ///
    /// `yindex` : `&usize`
    /// - index of the y location
    ///
    /// # Returns
    /// `Option<Vec<(usize, usize)>>`
    /// - `Some(Vec<(usize, usize)>)` : corner points in surrounding given
    ///   `x_index` and `y_index` in clockwise order.
    /// - `None` : `x_index` or `y_index` is out of range and no value exists.
    fn four_corners(&self, xindex: &usize, yindex: &usize) -> Option<Vec<(usize, usize)>> {
        if *xindex == 0
            || *yindex == 0
            || *xindex >= self.x.len() - 1
            || *yindex >= self.y.len() - 1
        {
            return None;
        }
        // clockwise in order
        Some(vec![
            (*xindex, yindex + 1),
            (xindex + 1, *yindex),
            (*xindex, yindex - 1),
            (xindex - 1, *yindex),
        ])
    }

    /// Interpolate the depth using crate::interpolator::bilinear
    ///
    /// First, the index points are converted to the x and y values at those
    /// indexes, then the depth at that index is taken. Finally, these are used
    /// as arguments to `interpolator::bilinear`.
    ///
    /// # Arguments
    /// `index_points`: `&Vec<(usize, usize)>`
    /// - a vector of (x_index, y_index) points representing the indices of the
    ///   corners that the target location is within.
    ///
    /// `target`: `&(f32, f32)`
    /// - interpolate the depth at this (x, y) point
    ///
    /// # Returns
    /// `Result<f32, Error>`
    /// - `Ok(f32)` : the depth at the target point
    /// - `Err(Error)` : cannot read depths from at coordinates in the `points`
    ///   vector.
    ///
    /// # Errors
    /// - `Error::IndexOutOfBounds` : one or more of the points passed to
    /// `points` is out of bounds.
    /// - `Error::InvalidArgument` : error during execution of
    /// `interpolator::bilinear` due to invalid arguments.
    fn interpolate(
        &self,
        index_points: &[(usize, usize)],
        target_point: &(f32, f32),
    ) -> Result<f32, Error> {
        let depth_points = vec![
            (
                self.x[index_points[0].0],
                self.y[index_points[0].1],
                self.depth_at_indexes(&index_points[0].0, &index_points[0].1)?,
            ),
            (
                self.x[index_points[1].0],
                self.y[index_points[1].1],
                self.depth_at_indexes(&index_points[1].0, &index_points[1].1)?,
            ),
            (
                self.x[index_points[2].0],
                self.y[index_points[2].1],
                self.depth_at_indexes(&index_points[2].0, &index_points[2].1)?,
            ),
            (
                self.x[index_points[3].0],
                self.y[index_points[3].1],
                self.depth_at_indexes(&index_points[3].0, &index_points[3].1)?,
            ),
        ];
        interpolator::bilinear(&depth_points, target_point)
    }

    /// Access values in flattened array as you would a 2d array
    ///
    /// # Arguments
    /// `x_index` : `&usize`
    /// - index of location in x array (column)
    ///
    /// `y_index` : `&usize`
    /// - index of location in y array (row)
    ///
    /// # Returns
    /// `Result<f32, Error>`
    /// - `Ok(f32)` : depth
    /// - `Err(Error::IndexOutOfBounds)` : the combined index (x_length *
    ///   y_index + x_index) is out of bounds of the depth array.
    ///
    /// # Errors
    /// `Err(Error::IndexOutOfBounds)` : this error is returned when `x_index`
    /// and `y_index` produce a value outside of the depth array.
    fn depth_at_indexes(&self, xindex: &usize, yindex: &usize) -> Result<f32, Error> {
        let index = self.x.len() * yindex + xindex;
        if index >= self.depth.len() {
            return Err(Error::IndexOutOfBounds);
        }
        Ok(self.depth[index])
    }
}

#[cfg(test)]
mod test_cartesian_file {

    use tempfile::NamedTempFile;

    use crate::{
        bathymetry::{cartesian_netcdf3::CartesianNetcdf3, BathymetryData},
        error::Error,
        io::utility::create_netcdf3_bathymetry,
    };

    /// create a file with four quadrants each with a different depth
    fn four_depth_fn(x: f32, y: f32) -> f64 {
        if x < 25000.0 {
            if y < 12500.0 {
                20.0
            } else {
                10.0
            }
        } else {
            if y < 12500.0 {
                5.0
            } else {
                15.0
            }
        }
    }

    #[test]
    // test accessing and viewing variables
    fn test_vars() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_bathymetry(&temp_path, 101, 51, 500.0, 500.0, four_depth_fn);

        let data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth");
        assert!((data.x[10] - 5000.0).abs() < f32::EPSILON)
    }

    #[test]
    // test the and view the nearest function
    fn test_get_nearest() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_bathymetry(&temp_path, 101, 51, 500.0, 500.0, four_depth_fn);

        let data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth");
        assert!(data.nearest(&5499.0, &data.x) == 11);
    }

    #[test]
    // check the output from four_corners function
    fn test_get_corners() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_bathymetry(&temp_path, 101, 51, 500.0, 500.0, four_depth_fn);

        let data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth");
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
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_bathymetry(&temp_path, 101, 51, 500.0, 500.0, four_depth_fn);

        let data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth");

        // check to see if depth is the same as above
        let check_depth = vec![
            (10099.0, 5099.0, 20.0),
            (30099.0, 5099.0, 5.0),
            (10099.0, 15099.0, 10.0),
            (30099.0, 15099.0, 15.0),
        ];

        for (x, y, h) in &check_depth {
            let depth = data.depth_and_gradient(x, y).unwrap().0;
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
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_bathymetry(&temp_path, 101, 51, 500.0, 500.0, four_depth_fn);

        let data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth");
        if let Error::IndexOutOfBounds = data.depth(&-500.1, &500.1).unwrap_err() {
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
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_bathymetry(&temp_path, 101, 51, 500.0, 500.0, four_depth_fn);

        let data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth");
        if let Error::IndexOutOfBounds = data.depth(&500.1, &-500.1).unwrap_err() {
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
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_bathymetry(&temp_path, 101, 51, 500.0, 500.0, four_depth_fn);

        let data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth");

        // check to see if depth is the same as above
        let check_depth = vec![
            (23000.0, 20000.0, 10.0),
            (10000.0, 12500.0, 12.5),
            (25000.0, 5000.0, 8.75),
            (40000.0, 12500.0, 12.5),
        ];

        for (x, y, h) in &check_depth {
            let depth = data.depth_and_gradient(x, y).unwrap().0;
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
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_bathymetry(&temp_path, 101, 51, 500.0, 500.0, four_depth_fn);

        let data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth");

        let nan = f32::NAN;

        assert!(data.depth(&nan, &nan).unwrap().is_nan());
        assert!(data.depth(&10000.0, &nan).unwrap().is_nan());
        assert!(data.depth(&nan, &10000.0).unwrap().is_nan());
    }

    #[test]
    fn test_depth_and_gradient_x() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        fn depth_fn(x: f32, _y: f32) -> f64 {
            x as f64 * 0.05
        }

        create_netcdf3_bathymetry(&temp_path, 100, 100, 1.0, 1.0, depth_fn);

        let data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth");

        // check to see if depth is the same as above
        let check_depth = vec![(10.0, 30.0, 0.5), (30.0, 10.0, 1.5)];

        for (x, y, h) in &check_depth {
            let depth = data.depth_and_gradient(x, y).unwrap().0;
            assert!(
                (depth - h).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                h,
                depth
            );
        }

        // check to see if gradient is the same
        let check_gradient = vec![
            (50.0, 50.0, 0.05, 0.0),
            (14.0, 12.0, 0.05, 0.0),
            (10.0, 80.0, 0.05, 0.0),
        ];

        for (x, y, dhdx, dhdy) in &check_gradient {
            let x_gradient = data.depth_and_gradient(x, y).unwrap().1 .0;
            let y_gradient = data.depth_and_gradient(x, y).unwrap().1 .1;
            assert!(
                (x_gradient - dhdx).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                dhdx,
                x_gradient
            );
            assert!(
                (y_gradient - dhdy).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                dhdy,
                y_gradient
            );
        }
    }

    #[test]
    fn test_depth_and_gradient_y() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        fn depth_fn(_x: f32, y: f32) -> f64 {
            y as f64 * 0.05
        }

        create_netcdf3_bathymetry(&temp_path, 100, 100, 1.0, 1.0, depth_fn);

        let data = CartesianNetcdf3::open(&temp_path, "x", "y", "depth");

        // check to see if depth is the same as above
        let check_depth = vec![(10.0, 30.0, 1.5), (30.0, 10.0, 0.5)];

        for (x, y, h) in &check_depth {
            let depth = data.depth_and_gradient(x, y).unwrap().0;
            assert!(
                (depth - h).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                h,
                depth
            );
        }

        // check to see if gradient is the same
        let check_gradient = vec![
            (50.0, 50.0, 0.0, 0.05),
            (14.0, 12.0, 0.0, 0.05),
            (10.0, 80.0, 0.0, 0.05),
        ];

        for (x, y, dhdx, dhdy) in &check_gradient {
            let x_gradient = data.depth_and_gradient(x, y).unwrap().1 .0;
            let y_gradient = data.depth_and_gradient(x, y).unwrap().1 .1;
            assert!(
                (x_gradient - dhdx).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                dhdx,
                x_gradient
            );
            assert!(
                (y_gradient - dhdy).abs() < f32::EPSILON,
                "Expected {}, but got {}",
                dhdy,
                y_gradient
            );
        }
    }
}
