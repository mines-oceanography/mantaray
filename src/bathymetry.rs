//! Bathymetry

use crate::error::Error;

/// A trait used to give the function get_depth
pub(crate) trait BathymetryData {
    /// Returns the nearest depth for the given x, y coordinate.
    fn get_depth(&self, x: &f32, y: &f32) -> Result<f32, Error>;
}

/// Read data from test_bathy_3.nc netcdf3 file that contains x, y, and depth
mod cartesian {

    use std::path::Path;
    use netcdf3::FileReader;
    use crate::{interpolator, bathymetry::BathymetryData, error::Error};

    /// A struct that stores a netcdf3 dataset named test_bathy_3.nc with
    /// methods to access, find nearest values, interpolate, and return depth.
    /// 
    /// # Note
    /// Currently, the methods do not know the difference between an out of
    /// bounds point and a just inside the bounds point. The nearest to each of
    /// these will be on the edge, so it will return an out of bounds error for
    /// both. In the future, the methods should be able to distinguish these two
    /// cases.
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
        /// - `Error::IndexOutOfBounds` : this error is returned when the
        /// `x` and `y` input give an out of bounds output.
        /// - `Error::InvalidArgument` : Currently
        /// as hard-coded, this function should never return this error.
        /// 
        /// # Panics
        /// This should only panic if the inputted `x` and `y` are out of
        /// bounds.
        fn get_depth(&self, x: &f32, y: &f32) -> Result<f32, Error> {
            self.depth(x, y)
        }
    }

    impl CartesianFile {
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
                data.read_var_f32("y").unwrap() ,
                data.read_var_f64("depth").unwrap().iter().map(|x| *x as f32).collect()
            );
            CartesianFile { variables }
        }
        
        /// Find nearest point
        /// 
        /// # Arguments
        /// `target` : `f32`
        /// - the value to find
        /// 
        /// `direction` : `&str`
        /// - `"x"` or `"y"`
        /// 
        /// # Returns
        /// `Result<usize, Error>`
        /// - `Ok(usize)` - index of closest value
        /// - `Err(Error)` - invalid `&str` passed for
        ///   `direction`
        /// 
        /// # Errors
        /// `Error::InvalidArgument` : this error is returned when `direction`
        /// is anything other than `"x"` or `"y"`.
        fn nearest(&self, target: &f32, direction: &str) -> Result<usize, Error> {
            let arr = match direction {
                "x" => &self.variables.0,
                "y" => &self.variables.1,
                _ => return Err(Error::InvalidArgument),
            };

            let mut closest_index = 0;
            let mut closest_distance = (target - arr[0]).abs();

            for i in 1..arr.len() {
                let distance = (target - arr[i]).abs();

                if distance < closest_distance {
                    closest_index = i;
                    closest_distance = distance;
                }
            }
            Ok(closest_index)
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
        /// `Result<(usize, usize), Error>`
        /// - `Ok((usize, usize))` : the nearest point to the given `x`, `y`
        ///   input.
        /// - `Error` : the passed `x` and `y` arguments create an index that is
        ///   out of bounds.
        /// 
        /// # Errors
        /// - `Error::IndexOutOfBounds` : this error is returned when the
        /// nearest_point is on the edge.
        /// - `Error::InvalidArgument` : Currently
        /// as hard-coded, this function should never return this error.
        /// 
        /// # Note
        /// This function will never panic, but if given an out of bounds point,
        /// it will return the closest edge. To attempt to fix this problem,
        /// `nearest_point` should error on points that are on the edges.
        fn nearest_point(&self, x: &f32, y: &f32) -> Result<(usize, usize), Error> {
            let indy = self.nearest(y, "y")?;
            let indx = self.nearest(x, "x")?;

            if indx <= 0 || indy <= 0 || indx >= self.variables.0.len()-1 || indy >= self.variables.1.len() {
                return Err(Error::IndexOutOfBounds);
            }

            Ok((indx, indy))
        }

        /// Get four adjecent points
        /// 
        /// # Arguments
        /// `indx` : `usize`
        /// - index of the x location
        /// 
        /// `indy` : `usize`
        /// - index of the y location
        /// 
        /// # Returns
        /// `Result<Vec<(usize, usize)>, Error>`
        /// - `Ok(Vec<(usize, usize)>)` : indices for the four corners
        /// surrounding the given indices in clockwise order.
        /// - `Err(Error::IndexOutOfBounds)` : the input indexes are out of
        ///   bounds of the depth array.
        /// 
        /// # Errors
        /// `Err(Error::IndexOutOfBounds)` : either `indx` or `indy` are out of bounds of the array.
        fn four_corners(&self, indx: &usize, indy: &usize) -> Result<Vec<(usize, usize)>, Error> {
            if *indx <= 0 || *indy <= 0 || 
            *indx >= self.variables.0.len()-1 || *indy >= self.variables.1.len()-1 {
                return Err(Error::IndexOutOfBounds);
            }
            // clockwise in order
            Ok(
                vec![
                (*indx, indy+1),
                (indx+1, *indy),
                (*indx, indy-1),
                (indx-1, *indy)
            ]
        )
        }
        
        /// Interpolate the depth using crate::interpolator::bilinear
        /// 
        /// # Arguments
        /// `points`: `&Vec<(usize, usize)>`
        /// - a vector of defined points in the depth grid
        /// 
        /// `target`: `&(f32, f32)`
        /// - interpolate this point
        /// 
        /// # Returns
        /// `Result<f32, Error>`
        /// - `Ok(f32)` : the depth at the target point
        /// - `Err(Error)` : cannot read depths from at coordinates in the
        ///   `points` vector.
        /// 
        /// # Errors
        /// `Error::IndexOutOfBounds` : one or more of the points passed to
        /// `points` is out of bounds.
        /// 
        /// # Panics
        /// This function will panic if `interpolator::bilinear` panics.
        fn interpolate(&self, points: &Vec<(usize, usize)>, target: &(f32, f32)) -> Result<f32, Error> {
            let pts = vec![
                (points[0].0 as i32, points[0].1 as i32, self.depth_from_arr(&points[0].0, &points[0].1)?),
                (points[1].0 as i32, points[1].1 as i32, self.depth_from_arr(&points[1].0, &points[1].1)?),
                (points[2].0 as i32, points[2].1 as i32, self.depth_from_arr(&points[2].0, &points[2].1)?),
                (points[3].0 as i32, points[3].1 as i32, self.depth_from_arr(&points[3].0, &points[3].1)?),
            ];
            Ok(interpolator::bilinear(&pts, target)?)
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
        
        /// Return the depth at x, y
        /// 
        /// Calculates the nearest point, finds surrounding edge points, then
        /// interpolates the depth at x,y.
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
        /// - `Err(Error)` : 
        /// 
        /// # Errors
        /// - `Error::IndexOutOfBounds` : this error is returned when the
        /// `x` and `y` input give an out of bounds output.
        /// - `Error::InvalidArgument` : Currently
        /// as hard-coded, this function should never return this error.
        /// 
        /// # Panics
        /// This will only panic if the `interpolator::bilinear` function in `self.interpolate` panics.
        fn depth(&self, x: &f32, y: &f32) -> Result<f32, Error> {
            let nearest_pt = self.nearest_point(x, y)?;
            let edge_points = self.four_corners(&nearest_pt.0, &nearest_pt.1)?;
            let depth = self.interpolate(&edge_points, &(*x, *y))?;
            Ok(depth)
        }

    }

    #[test]
    // test accessing and viewing variables
    fn test_vars() {
        let data = CartesianFile::new(Path::new("data/test_bathy_3.nc"));
        dbg!(data.variables.0[10]);
    }

    #[test]
    // test the and view the nearest function
    fn test_get_nearest() {
        let data = CartesianFile::new(Path::new("data/test_bathy_3.nc"));
        dbg!(data.nearest(&5499.0, "x").unwrap());
        assert!(data.nearest(&5499.0, "x").unwrap() == 11);
    }

    #[test]
    // view the output from four_corners function
    fn test_get_corners() {
        let data = CartesianFile::new(Path::new("data/test_bathy_3.nc"));
        dbg!(data.four_corners(&10, &10).unwrap());
    }

    #[test]
    // check values inside the four quadrants but not on grid points
    fn test_depth() {
        let data = CartesianFile::new(Path::new("data/test_bathy_3.nc"));
        assert!(data.get_depth(&10099.0, &5099.0).unwrap() == 20.0);
        assert!(data.get_depth(&30099.0, &5090.0).unwrap() == 5.0);
        assert!(data.get_depth(&10099.0, &15099.0).unwrap() == 10.0);
        assert!(data.get_depth(&30099.0, &15099.0).unwrap() == 15.0);
    }

    #[test]
    fn test_out_of_bounds() {
        let data = CartesianFile::new(Path::new("data/test_bathy_3.nc"));
        if let Error::IndexOutOfBounds = data.get_depth(&-500.1, &-500.1).unwrap_err() {
            assert!(true);
        } else {
            assert!(false);
        }

    }

}