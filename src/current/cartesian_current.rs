//! Handle netcdf files in cartesian coordinates containing snapshot of current
//! field.

use std::path::Path;

use netcdf3::{DataType, FileReader};

use super::CurrentData;
use crate::error::Error;
use crate::error::Result;
use crate::interpolator;
use crate::{Current, Point};

#[derive(Debug)]
#[allow(dead_code)]
/// A struct to hold the data from a NetCDF file in a Cartesian coordinates with
/// x, y, u, and v values constant in time.
pub(crate) struct CartesianCurrent {
    /// vector of the x variable
    x_vec: Vec<f64>,
    /// vector of the y variable
    y_vec: Vec<f64>,
    /// vector of the u variable
    u_vec: Vec<f64>,
    /// vector of the v variable
    v_vec: Vec<f64>,
}

#[allow(dead_code)]
impl CartesianCurrent {
    /// Create a new `Cartesian` from a NetCDF file.
    ///
    /// # Arguments
    /// - `path` : `&Path` Path to the NetCDF file.
    ///
    /// - `x_name` : `&str` Name of the variable in the NetCDF file that
    ///   contains the x data.
    ///
    /// - `y_name` : `&str` Name of the variable in the NetCDF file that
    ///   contains the y data.
    ///
    /// - `u_name` : `&str` Name of the variable in the NetCDF file that
    ///   contains the u data.
    ///
    /// - `v_name` : `&str` Name of the variable in the NetCDF file that
    ///   contains the v data.
    ///
    /// # Returns
    /// `Self` : `CurrentCartesianFile` the new constructed struct.
    ///
    /// # Panics
    /// Panics if the NetCDF file does not contain the variables `x`, `y`, `u`,
    /// `v`.
    ///
    /// # Note
    /// The variables `x`, `y`, `u`, `v` can be of any type that is in
    /// `netcdf3::DataType`.
    pub(crate) fn open(
        path: &Path,
        x_name: &str,
        y_name: &str,
        u_name: &str,
        v_name: &str,
    ) -> Self {
        let mut data = FileReader::open(path).unwrap();

        let x_data = data.read_var(x_name).unwrap();
        let x_data = match x_data.data_type() {
            DataType::I16 => x_data
                .get_i16_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::I8 => x_data
                .get_i8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::U8 => x_data
                .get_u8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::I32 => x_data
                .get_i32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::F32 => x_data
                .get_f32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::F64 => x_data.get_f64_into().unwrap(),
        };

        let y_data = data.read_var(y_name).unwrap();
        let y_data = match y_data.data_type() {
            DataType::I16 => y_data
                .get_i16_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::I8 => y_data
                .get_i8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::U8 => y_data
                .get_u8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::I32 => y_data
                .get_i32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::F32 => y_data
                .get_f32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::F64 => y_data.get_f64_into().unwrap(),
        };

        let u_data = data.read_var(u_name).unwrap();
        let u_data = match u_data.data_type() {
            DataType::I16 => u_data
                .get_i16_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::I8 => u_data
                .get_i8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::U8 => u_data
                .get_u8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::I32 => u_data
                .get_i32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::F32 => u_data
                .get_f32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::F64 => u_data.get_f64_into().unwrap(),
        };

        let v_data = data.read_var(v_name).unwrap();
        let v_data = match v_data.data_type() {
            DataType::I16 => v_data
                .get_i16_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::I8 => v_data
                .get_i8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::U8 => v_data
                .get_u8_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::I32 => v_data
                .get_i32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::F32 => v_data
                .get_f32_into()
                .unwrap()
                .iter()
                .map(|x| *x as f64)
                .collect(),
            DataType::F64 => v_data.get_f64_into().unwrap(),
        };

        CartesianCurrent {
            x_vec: x_data,
            y_vec: y_data,
            u_vec: u_data,
            v_vec: v_data,
        }
    }

    /// Find nearest point
    ///
    /// # Arguments
    /// `target` : `f64`
    /// - the value to find
    ///
    /// `arr` : `&[f64]`
    /// - the array that will be used when searching for the closest value.
    ///
    /// # Returns
    /// `usize`: index of closest value
    ///
    /// # Note
    /// Binary search might be faster, but this is a simple implementation. And
    /// binary search requires it to be sorted.
    fn nearest(&self, target: &f64, arr: &[f64]) -> usize {
        let mut left = 0;
        let mut right = arr.len() - 1;
        let mut closest_index = 0;
        let mut closest_distance = f64::MAX;

        // edge cases
        if *target <= arr[left] {
            return left;
        }
        if *target >= arr[right] {
            return right;
        }

        // binary search
        while left <= right {
            let mid = (left + right) / 2;
            let distance = (target - arr[mid]).abs();

            if distance < closest_distance {
                closest_index = mid;
                closest_distance = distance;
            }

            if arr[mid] == *target {
                return mid;
            }

            if arr[mid] > *target {
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }

        closest_index
    }

    /// Nearest (x, y) index point for given coordinate
    ///
    /// # Arguments
    /// `x`: `&f64`
    /// - x coordinate
    ///
    /// `y`: `&f64`
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
    fn nearest_point(&self, x: &f64, y: &f64) -> Option<(usize, usize)> {
        let indx = self.nearest(x, &self.x_vec);
        let indy = self.nearest(y, &self.y_vec);

        if indx == 0 || indy == 0 || indx >= self.x_vec.len() - 1 || indy >= self.y_vec.len() - 1 {
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
    ///
    /// NOTE: with the addition of the time dimension, this function will need
    /// to be updated to include the time dimension. Therefore, it will need to
    /// return a vec of 6 (t, x, y)
    fn four_corners(&self, indx: &usize, indy: &usize) -> Option<Vec<(usize, usize)>> {
        if *indx == 0
            || *indy == 0
            || *indx >= self.x_vec.len() - 1
            || *indy >= self.y_vec.len() - 1
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
    /// `points`: `&[(usize, usize)]`
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
    fn interpolate(
        &self,
        points: &[(usize, usize)], // 4 points
        target: &(f32, f32),
        value_arr: &[f64],
    ) -> Result<f32> {
        if points.len() != 4 {
            return Err(Error::InvalidArgument);
        }

        let pts = vec![
            (
                self.x_vec[points[0].0] as f32,                                   // x1
                self.y_vec[points[0].1] as f32,                                   // y1
                self.val_from_arr(&points[0].0, &points[0].1, value_arr)? as f32, // z1
            ),
            (
                self.x_vec[points[1].0] as f32,
                self.y_vec[points[1].1] as f32,
                self.val_from_arr(&points[1].0, &points[1].1, value_arr)? as f32,
            ),
            (
                self.x_vec[points[2].0] as f32,
                self.y_vec[points[2].1] as f32,
                self.val_from_arr(&points[2].0, &points[2].1, value_arr)? as f32,
            ),
            (
                self.x_vec[points[3].0] as f32,
                self.y_vec[points[3].1] as f32,
                self.val_from_arr(&points[3].0, &points[3].1, value_arr)? as f32,
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
    /// `arr` : `&[f64]`
    /// - the array to access
    ///
    /// # Returns
    /// `Result<f64, Error>`
    /// - `Ok(f64)` : value at the given index
    /// - `Err(Error::IndexOutOfBounds)` : the combined index (x_length *
    ///   indy + indx) is out of bounds of array.
    ///
    /// # Errors
    /// `Err(Error::IndexOutOfBounds)` : this error is returned when `indx`
    /// and `indy` produce a value outside of the array.
    fn val_from_arr(&self, indx: &usize, indy: &usize, arr: &[f64]) -> Result<f64> {
        let index = self.x_vec.len() * indy + indx;
        if index >= arr.len() {
            return Err(Error::IndexOutOfBounds);
        }
        Ok(arr[index])
    }
}

impl CurrentData for CartesianCurrent {
    /// return the current at the point (x, y)
    ///
    /// # Arguments
    ///
    /// - `x` : `&f64` x coordinate
    ///
    /// - `y` : `&f64` y coordinate
    ///
    /// # Returns
    ///
    /// `Result<(f64, f64), Error>` : the current at the point (x, y) or an
    /// error
    ///
    /// # Errors
    ///
    /// `Error::IndexOutOfBounds` : the point (x, y) is out of bounds of the
    /// data
    fn current(&self, point: &Point<f64>) -> Result<Current<f64>> {
        // get the nearest point
        let (indx, indy) = match self.nearest_point(point.x(), point.y()) {
            Some((indx, indy)) => (indx, indy),
            None => return Err(Error::IndexOutOfBounds),
        };

        // get the four corners
        let corners = match self.four_corners(&indx, &indy) {
            Some(corners) => corners,
            None => return Err(Error::IndexOutOfBounds),
        };

        // interpolate the u and v values
        let u = self.interpolate(
            &corners,
            &(*point.x() as f32, *point.y() as f32),
            &self.u_vec,
        )?;
        let v = self.interpolate(
            &corners,
            &(*point.x() as f32, *point.y() as f32),
            &self.v_vec,
        )?;

        Ok(Current::new(u as f64, v as f64))
    }

    /// return the current and the gradient at the point (x, y)
    ///
    /// # Arguments
    ///
    /// - `x` : `&f64` x coordinate
    ///
    /// - `y` : `&f64` y coordinate
    ///
    /// # Returns
    ///
    /// `Result<((f64, f64), (f64, f64, f64, f64)), Error>` : the current at the
    /// point (x, y) and the gradient at the point (x, y) or an error.
    ///
    /// # Errors
    ///
    /// `Error::IndexOutOfBounds` : the point (x, y) is out of bounds of the
    /// data
    fn current_and_gradient(
        &self,
        point: &Point<f64>,
    ) -> Result<((f64, f64), (f64, f64, f64, f64))> {
        // get the nearest point
        let (indx, indy) = match self.nearest_point(point.x(), point.y()) {
            Some((indx, indy)) => (indx, indy),
            None => return Err(Error::IndexOutOfBounds),
        };

        // get the four corners
        let corners = match self.four_corners(&indx, &indy) {
            Some(corners) => corners,
            None => return Err(Error::IndexOutOfBounds),
        };

        // interpolate the u and v values
        let u = self.interpolate(
            &corners,
            &(*point.x() as f32, *point.y() as f32),
            &self.u_vec,
        )?;
        let v = self.interpolate(
            &corners,
            &(*point.x() as f32, *point.y() as f32),
            &self.v_vec,
        )?;

        // calculate the gradients

        // NOTE: the gradient assumes that the depth is linear in both the x
        // and y directions, and since bilinear interpolation is used to
        // interpolate the depth at any given point, this is a good
        // approximation.
        let x_space = self.x_vec[1] - self.x_vec[0];
        let y_space = self.y_vec[1] - self.y_vec[0];

        let dudx = (self.val_from_arr(&corners[1].0, &corners[1].1, &self.u_vec)?
            - self.val_from_arr(&corners[3].0, &corners[3].1, &self.u_vec)?)
            / (2.0 * x_space);

        let dudy = (self.val_from_arr(&corners[0].0, &corners[0].1, &self.u_vec)?
            - self.val_from_arr(&corners[2].0, &corners[2].1, &self.u_vec)?)
            / (2.0 * y_space);

        let dvdx = (self.val_from_arr(&corners[1].0, &corners[1].1, &self.v_vec)?
            - self.val_from_arr(&corners[3].0, &corners[3].1, &self.v_vec)?)
            / (2.0 * x_space);

        let dvdy = (self.val_from_arr(&corners[0].0, &corners[0].1, &self.v_vec)?
            - self.val_from_arr(&corners[2].0, &corners[2].1, &self.v_vec)?)
            / (2.0 * y_space);

        Ok(((u as f64, v as f64), (dudx, dudy, dvdx, dvdy)))
    }
}

#[cfg(test)]
mod test_cartesian_file_current {
    use tempfile::NamedTempFile;

    use crate::{
        current::{cartesian_current::CartesianCurrent, CurrentData},
        io::utility::create_netcdf3_current,
        Current, Point,
    };
    use std::path::Path;

    /// returns a simple current with u = 5 and v = 0
    fn simple_current(_x: f32, _y: f32) -> (f64, f64) {
        (5.0, 0.0)
    }

    /// this will create a current file it will have x and y as f32 and u and v
    /// as f64. this will have a gradient in the u and v fields
    fn simple_x_gradient(x: f32, _y: f32) -> (f64, f64) {
        let x = x as f64;
        (x, x)
    }

    /// this will create a current file it will have x and y as f32 and u and v
    /// as f64. this will have a gradient in the u and v fields
    fn simple_y_gradient(_x: f32, y: f32) -> (f64, f64) {
        let y = y as f64;
        (y, y)
    }

    /// create a current file with variable (x, y) as (i16, i8) and (u, v) as
    /// (u8, i32). this is a special case file just for testing purposes, so it
    /// stays for now.
    fn create_netcdf3_current_iu(
        path: &Path,
        x_len: usize,
        y_len: usize,
        x_step: f32,
        y_step: f32,
    ) {
        let x_data: Vec<i16> = (0..x_len).map(|x| x as i16 * x_step as i16).collect();
        let y_data: Vec<i8> = (0..y_len).map(|y| y as i8 * y_step as i8).collect();

        let u_data: Vec<u8> = (0..x_len * y_len).map(|_| 5_u8).collect();
        let v_data: Vec<i32> = (0..x_len * y_len).map(|_| 0_i32).collect();

        // most below copied from the docs
        use netcdf3::{DataSet, FileWriter, Version};
        let y_dim_name: &str = "y";
        let y_var_name: &str = y_dim_name;
        let y_var_len: usize = y_len;

        let x_dim_name: &str = "x";
        let x_var_name: &str = x_dim_name;
        let x_var_len: usize = x_len;

        let u_dim_name: &str = "u";
        let u_var_name: &str = u_dim_name;
        let u_var_len: usize = u_data.len();

        let v_dim_name: &str = "v";
        let v_var_name: &str = v_dim_name;
        let v_var_len: usize = v_data.len();

        // Create the NetCDF-3 definition
        // ------------------------------
        assert_eq!(u_var_len, y_var_len * x_var_len);
        assert_eq!(v_var_len, y_var_len * x_var_len);
        let data_set: DataSet = {
            let mut data_set: DataSet = DataSet::new();
            // Define the dimensions
            data_set.add_fixed_dim(y_dim_name, y_var_len).unwrap();
            data_set.add_fixed_dim(x_dim_name, x_var_len).unwrap();
            // Define the variable
            data_set.add_var_i8(y_var_name, &[y_dim_name]).unwrap();
            data_set.add_var_i16(x_var_name, &[x_var_name]).unwrap();
            data_set
                .add_var_u8(u_var_name, &[y_dim_name, x_var_name])
                .unwrap();
            data_set
                .add_var_i32(v_var_name, &[y_dim_name, x_var_name])
                .unwrap();

            data_set
        };

        // Create and write the NetCDF-3 file
        // ----------------------------------
        let mut file_writer: FileWriter = FileWriter::open(path).unwrap();
        // Set the NetCDF-3 definition
        file_writer.set_def(&data_set, Version::Classic, 0).unwrap();
        file_writer.write_var_i8(y_var_name, &y_data[..]).unwrap();
        file_writer.write_var_i16(x_var_name, &x_data[..]).unwrap();
        file_writer.write_var_u8(u_var_name, &u_data[..]).unwrap();
        file_writer.write_var_i32(v_var_name, &v_data[..]).unwrap();
        // file_writer.close().unwrap();
        // end of copied from docs
    }

    #[test]
    fn test_all_types() {
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        // test with f32 and f64
        create_netcdf3_current(&path, 1, 1, 1.0, 1.0, simple_current);
        let _: CartesianCurrent = CartesianCurrent::open(&path, "x", "y", "u", "v");

        // test with i16, i8, u8, i32
        create_netcdf3_current_iu(&path, 1, 1, 1.0, 1.0);
        let _: CartesianCurrent = CartesianCurrent::open(&path, "x", "y", "u", "v");
    }

    #[test]
    // test the and view the nearest function. Note that this is the same test
    // case as in bathymetry/cartesian.rs
    fn test_get_nearest() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 100, 100, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");

        // inside the bounds
        assert_eq!(data.nearest(&5499.0, &data.x_vec), 11);

        // test out of bounds
        assert_eq!(data.nearest(&-5499.0, &data.x_vec), 0);
        assert_eq!(data.nearest(&50_001.0, &data.x_vec), 99);
    }

    #[test]
    // test the and nearest point function on one point
    fn test_get_nearest_point() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 101, 51, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");

        // inside the bounds
        assert!(data.nearest_point(&5499.0, &499.0) == Some((11, 1)));

        // test out of bounds
        assert!(data.nearest_point(&-5499.0, &-499.0) == None);
        assert!(data.nearest_point(&-5499.0, &50_001.0) == None);
        assert!(data.nearest_point(&50_001.0, &50_001.0) == None);
    }

    #[test]
    // check the output from four_corners function
    fn test_get_corners() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 101, 51, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");

        // in bounds
        let corners = data.four_corners(&10, &10).unwrap();
        assert!(corners[0].0 == 10 && corners[0].1 == 11);
        assert!(corners[1].0 == 11 && corners[1].1 == 10);
        assert!(corners[2].0 == 10 && corners[2].1 == 9);
        assert!(corners[3].0 == 9 && corners[3].1 == 10);

        // out of bounds
        let corners = data.four_corners(&0, &0);
        assert!(corners == None);

        let corners = data.four_corners(&100, &100);
        assert!(corners == None);
    }

    #[test]
    // test the interpolate function
    fn test_interpolate() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 101, 51, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");
        let corners = data.four_corners(&10, &10).unwrap();
        let interpolated = data.interpolate(&corners, &(5499.0, 499.0), &data.u_vec);
        assert!(interpolated.unwrap() == 5.0);

        let interpolated = data.interpolate(&corners, &(5499.0, 499.0), &data.v_vec);
        assert!(interpolated.unwrap() == 0.0);
    }

    #[test]
    // test the value_from_arr function
    fn test_val_from_arr() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 101, 51, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");
        let val = data.val_from_arr(&10, &10, &data.u_vec);
        assert!(val.unwrap() == 5.0);

        let val = data.val_from_arr(&10, &10, &data.v_vec);
        assert!(val.unwrap() == 0.0);

        // test out of bounds
        let val = data.val_from_arr(&100, &100, &data.u_vec);
        assert!(val.is_err());
    }

    #[test]
    // test the current function
    fn test_current() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 101, 51, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");
        let current = data.current(&Point::new(5499.0, 499.0));
        assert!(current.unwrap() == Current::new(5.0, 0.0));

        // test out of bounds
        let current = data.current(&Point::new(50_001.0, 1000.0));
        assert!(current.is_err());

        let current = data.current(&Point::new(-50_001.0, -1000.0));
        assert!(current.is_err());
    }

    #[test]
    // test the current_and_gradient function
    fn test_current_and_zero_grad() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 101, 51, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");
        let current = data.current_and_gradient(&Point::new(5499.0, 499.0));
        assert!(current.unwrap() == ((5.0, 0.0), (0.0, 0.0, 0.0, 0.0)));

        // test out of bounds
        let current = data.current_and_gradient(&Point::new(50_001.0, 1000.0));
        assert!(current.is_err());

        let current = data.current_and_gradient(&Point::new(-50_001.0, -1000.0));
        assert!(current.is_err());
    }

    #[test]
    // test the current_and_gradient function with constant gradients in x direction
    fn test_current_and_grad_x() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 100, 100, 1.0, 1.0, simple_x_gradient);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");
        let current = data.current_and_gradient(&Point::new(45.0, 45.0));
        assert_eq!(current.unwrap(), ((45.0, 45.0), (1.0, 0.0, 1.0, 0.0)));
    }

    #[test]
    // test the current_and_gradient function with constant gradients in y direction
    fn test_current_and_grad_y() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 100, 100, 1.0, 1.0, simple_y_gradient);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");
        let current = data.current_and_gradient(&Point::new(45.0, 45.0));
        assert_eq!(current.unwrap(), ((45.0, 45.0), (0.0, 1.0, 0.0, 1.0)));
    }
}
