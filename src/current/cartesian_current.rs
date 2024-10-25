//! Handle netcdf files in cartesian coordinates containing snapshot of current
//! field.

use std::path::Path;

use netcdf3::{DataType, FileReader};

use super::CurrentData;
use crate::datatype::{Current, Gradient, Point};
use crate::error::Error;
use crate::error::Result;
use crate::interpolator;

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

    /// Find the index of the closest value to the target in the array
    ///
    /// # Arguments
    /// `target` : `&f64`
    /// - the value to find
    ///
    /// `arr` : `&[f64]`
    /// - the array that will be used when searching for the closest value.
    ///
    /// # Returns
    /// `Result<f64>`: index of closest value or error
    ///
    /// # Note
    /// This function assumes the array has equal spacing between all elements
    /// and is ordered from least to greatest. Given those two conditions, it is
    /// valid to have fractional indexes.
    fn nearest(&self, target: &f64, array: &[f64]) -> Result<f64> {
        // array has to have at least 1 element (prevent future divide by zero error)
        if array.is_empty() {
            return Err(Error::IndexOutOfBounds); // error
        }

        // if the array has only one element, return 0 as its the only option
        if array.len() == 1 {
            return Ok(0.0);
        }

        // we know the array has at least two elements, so the following line
        // will never panic
        let spacing = (array[1] - array[0]).abs();

        let index = (target - array[0]) / spacing;

        if index < 0.0 || index > (array.len() - 1) as f64 {
            Err(Error::IndexOutOfBounds)
        } else {
            Ok(index)
        }
    }

    /// Returns the nearest (xindex, yindex) point to given (x ,y) point
    ///
    /// # Arguments
    /// `point` : `&Point<f64>` the location of the ray
    ///
    /// # Returns
    /// `Result<(f64, f64)>`: the indexes of the nearest point or an error.
    ///
    /// # Note
    /// This function assumes the x and y dimensions of the data are equally
    /// spaced arrays in ascending order. Therefore, fractional indexes are expected.
    fn nearest_point(&self, point: &Point<f64>) -> Result<(f64, f64)> {
        // find floating point "index"
        let xindex = self.nearest(point.x(), &self.x_vec)?;
        let yindex = self.nearest(point.y(), &self.y_vec)?;

        Ok((xindex, yindex))
    }

    /// Get four adjacent points
    ///
    /// # Arguments
    /// `point` : `&Point<f64>` the point to find the 4 corners around
    ///
    /// # Returns
    /// `Result<Vec<(usize, usize)>>`: returns a vector of the 4 points
    /// surrounding the target point. The points are in clockwise order starting
    /// with the bottom left point. Or it will return an out of bounds error.
    fn four_corners(&self, point: &Point<f64>) -> Result<Vec<(usize, usize)>> {
        let (xindex, yindex) = self.nearest_point(point)?;

        // determine the edges
        let xlow = 0.0;
        let xhigh = (self.x_vec.len() - 1) as f64;
        let ylow = 0.0;
        let yhigh = (self.y_vec.len() - 1) as f64;

        // check edges, interior points, or normal case
        let (x1, x2) = if xindex == xlow {
            // left edge
            let x1 = xindex as usize;
            let x2 = xindex as usize + 1;
            (x1, x2)
        } else if xindex == xhigh {
            // right edge
            let x1 = xindex as usize - 1;
            let x2 = xindex as usize;
            (x1, x2)
        } else if xindex.fract() == 0.0 {
            // on x grid point, but not on edge
            let x1 = xindex.round() as usize;
            let x2 = x1 + 1;
            (x1, x2)
        } else {
            // normal case
            let x1 = xindex.floor() as usize;
            let x2 = xindex.ceil() as usize;
            (x1, x2)
        };

        // check edges, interior points, or normal case
        let (y1, y2) = if yindex == ylow {
            // bottom edge
            let y1 = yindex as usize;
            let y2 = yindex as usize + 1;
            (y1, y2)
        } else if yindex == yhigh {
            // top edge
            let y1 = yindex as usize - 1;
            let y2 = yindex as usize;
            (y1, y2)
        } else if yindex.fract() == 0.0 {
            // on y grid point, but not edge
            let y1 = yindex.round() as usize;
            let y2 = y1 + 1;
            (y1, y2)
        } else {
            // normal case
            let y1 = yindex.floor() as usize;
            let y2 = yindex.ceil() as usize;
            (y1, y2)
        };

        Ok(vec![(x1, y1), (x1, y2), (x2, y2), (x2, y1)])
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
        // get the four corners
        let corners = match self.four_corners(point) {
            Ok(corners) => corners,
            Err(e) => return Err(e),
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
    ) -> Result<(Current<f64>, (Gradient<f64>, Gradient<f64>))> {
        // get the four corners
        let corners = match self.four_corners(point) {
            Ok(corners) => corners,
            Err(e) => return Err(e),
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

        let sw_point = &corners[0];
        let nw_point = &corners[1];
        let se_point = &corners[3];

        let dudx = (self.val_from_arr(&se_point.0, &se_point.1, &self.u_vec)?
            - self.val_from_arr(&sw_point.0, &sw_point.1, &self.u_vec)?)
            / x_space;

        let dudy = (self.val_from_arr(&nw_point.0, &nw_point.1, &self.u_vec)?
            - self.val_from_arr(&sw_point.0, &sw_point.1, &self.u_vec)?)
            / y_space;

        let dvdx = (self.val_from_arr(&se_point.0, &se_point.1, &self.v_vec)?
            - self.val_from_arr(&sw_point.0, &sw_point.1, &self.v_vec)?)
            / x_space;

        let dvdy = (self.val_from_arr(&nw_point.0, &nw_point.1, &self.v_vec)?
            - self.val_from_arr(&sw_point.0, &sw_point.1, &self.v_vec)?)
            / y_space;

        Ok((
            Current::new(u as f64, v as f64),
            (Gradient::new(dudx, dudy), Gradient::new(dvdx, dvdy)),
        ))
    }
}

#[cfg(test)]
mod test_cartesian_file_current {
    use tempfile::NamedTempFile;

    use super::{Current, Gradient, Point};
    use crate::{
        current::{cartesian_current::CartesianCurrent, CurrentData},
        error::Error,
        io::utility::create_netcdf3_current,
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
    // test the and view the nearest function
    fn test_nearest() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_current(&temp_path, 101, 51, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(&temp_path, "x", "y", "u", "v");

        // in bounds
        assert!(data.nearest(&5499.0, &data.x_vec).unwrap().round() == 11.0);

        // out of bounds
        assert!(data.nearest(&-1.0, &data.y_vec).is_err());
        assert!(data.nearest(&25_501.0, &data.y_vec).is_err());

        // on grid point
        assert!((data.nearest(&5500.0, &data.x_vec).unwrap() - 11.0).abs() <= f64::EPSILON);
    }

    #[test]
    // test the nearest point function (which returns floating point indexes)
    fn test_nearest_point() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_current(&temp_path, 101, 51, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(&temp_path, "x", "y", "u", "v");

        // in bounds
        assert!(
            data.nearest_point(&Point::new(1.0, 24_999.0))
                .unwrap()
                .0
                .round()
                == 0.0
        );
        assert!(
            data.nearest_point(&Point::new(1.0, 24_999.0))
                .unwrap()
                .1
                .round()
                == 50.0
        );

        // out of bounds
        assert!(data.nearest_point(&Point::new(1.0, 25_001.0)).is_err());
        assert!(data.nearest_point(&Point::new(-1.0, 25_000.0)).is_err());

        // grid points
        assert!(
            (data.nearest_point(&Point::new(0.0, 25_000.0)).unwrap().0 - 0.0).abs() <= f64::EPSILON
        );
        assert!(
            (data.nearest_point(&Point::new(0.0, 25_000.0)).unwrap().1 - 50.0).abs()
                <= f64::EPSILON
        );
    }

    #[test]
    // check all the cases for the output from the four_corners function
    fn test_get_corners() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = temp_file.into_temp_path();

        create_netcdf3_current(&temp_path, 101, 51, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(&temp_path, "x", "y", "u", "v");

        // check edge cases

        // top left corner
        assert!(
            data.four_corners(&Point::new(0.0, 25_000.0)).unwrap()
                == vec![(0, 49), (0, 50), (1, 50), (1, 49)]
        );

        // left edge
        assert!(
            data.four_corners(&Point::new(0.0, 5_500.0)).unwrap()
                == vec![(0, 11), (0, 12), (1, 12), (1, 11)]
        );

        // bottom left corner
        assert!(
            data.four_corners(&Point::new(0.0, 0.0)).unwrap()
                == vec![(0, 0), (0, 1), (1, 1), (1, 0)]
        );

        // top edge
        assert!(
            data.four_corners(&Point::new(5_500.0, 25_000.0)).unwrap()
                == vec![(11, 49), (11, 50), (12, 50), (12, 49)]
        );

        // bottom edge
        assert!(
            data.four_corners(&Point::new(5_500.0, 0.0)).unwrap()
                == vec![(11, 0), (11, 1), (12, 1), (12, 0)]
        );

        // top right corner
        assert!(
            data.four_corners(&Point::new(50_000.0, 25_000.0)).unwrap()
                == vec![(99, 49), (99, 50), (100, 50), (100, 49)]
        );

        // right edge
        assert!(
            data.four_corners(&Point::new(50_000.0, 5_500.0)).unwrap()
                == vec![(99, 11), (99, 12), (100, 12), (100, 11)]
        );

        // bottom right corner
        assert!(
            data.four_corners(&Point::new(50_000.0, 0.0)).unwrap()
                == vec![(99, 0), (99, 1), (100, 1), (100, 0)]
        );

        // check out of bounds
        assert!(match data.four_corners(&Point::new(50_001.0, 0.0)) {
            Err(Error::IndexOutOfBounds) => true,
            _ => false,
        });
        assert!(match data.four_corners(&Point::new(50_000.0, 25_001.0)) {
            Err(Error::IndexOutOfBounds) => true,
            _ => false,
        });
        assert!(match data.four_corners(&Point::new(-1.0, 0.0)) {
            Err(Error::IndexOutOfBounds) => true,
            _ => false,
        });
        assert!(match data.four_corners(&Point::new(50_000.0, -1.0)) {
            Err(Error::IndexOutOfBounds) => true,
            _ => false,
        });

        // check not edge, in bounds, and both x and y on grid point
        assert!(
            data.four_corners(&Point::new(5_500.0, 5_500.0)).unwrap()
                == vec![(11, 11), (11, 12), (12, 12), (12, 11)]
        );

        // check not edge, in bounds, and only x on grid point
        assert!(
            data.four_corners(&Point::new(5_500.0, 5_750.0)).unwrap()
                == vec![(11, 11), (11, 12), (12, 12), (12, 11)]
        );

        // check not edge, in bounds, and only y on grid point
        assert!(
            data.four_corners(&Point::new(5_750.0, 5_500.0)).unwrap()
                == vec![(11, 11), (11, 12), (12, 12), (12, 11)]
        );

        // check not edge, in bounds, and not on a grid point
        assert!(
            data.four_corners(&Point::new(5_750.0, 5_750.0)).unwrap()
                == vec![(11, 11), (11, 12), (12, 12), (12, 11)]
        );
    }

    #[test]
    // test the interpolate function
    fn test_interpolate() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 101, 51, 500.0, 500.0, simple_current);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");
        let corners = data.four_corners(&Point::new(10.0, 10.0)).unwrap();
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

        create_netcdf3_current(&path, 100, 50, 1.0, 1.0, simple_current);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");

        // check full domain is accurate
        for i in 0..100 {
            for j in 0..50 {
                let i = i as f64;
                let j = j as f64;

                let current = data.current(&Point::new(i, j)).unwrap();

                assert_eq!(current, Current::new(5.0, 0.0))
            }
        }

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

        create_netcdf3_current(&path, 100, 50, 1.0, 1.0, simple_current);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");

        // check full domain is accurate
        for i in 0..100 {
            for j in 0..50 {
                let i = i as f64;
                let j = j as f64;

                let current_and_gradient = data.current_and_gradient(&Point::new(i, j)).unwrap();

                assert_eq!(
                    current_and_gradient,
                    (
                        Current::new(5.0, 0.0),
                        (Gradient::new(0.0, 0.0), Gradient::new(0.0, 0.0))
                    )
                )
            }
        }

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

        // check full domain is accurate
        for i in 0..100 {
            for j in 0..100 {
                let i = i as f64;
                let j = j as f64;

                let current_and_gradient = data.current_and_gradient(&Point::new(i, j)).unwrap();

                assert_eq!(
                    current_and_gradient,
                    (
                        Current::new(i, i),
                        (Gradient::new(1.0, 0.0), Gradient::new(1.0, 0.0))
                    )
                )
            }
        }
    }

    #[test]
    // test the current_and_gradient function with constant gradients in y direction
    fn test_current_and_grad_y() {
        // create temporary file
        let temp_file = NamedTempFile::new().unwrap();
        let path = temp_file.into_temp_path();

        create_netcdf3_current(&path, 100, 100, 1.0, 1.0, simple_y_gradient);

        let data = CartesianCurrent::open(Path::new(&path), "x", "y", "u", "v");

        // check full domain is accurate
        for i in 0..100 {
            for j in 0..100 {
                let i = i as f64;
                let j = j as f64;

                let current_and_gradient = data.current_and_gradient(&Point::new(i, j)).unwrap();

                assert_eq!(
                    current_and_gradient,
                    (
                        Current::new(j, j),
                        (Gradient::new(0.0, 1.0), Gradient::new(0.0, 1.0))
                    )
                )
            }
        }
    }
}
