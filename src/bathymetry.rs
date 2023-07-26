//! Bathymetry

/// A trait used to give the function get_depth
pub(crate) trait BathymetryData {
    /// Returns the nearest depth for the given lat, lon coordinate.
    fn get_depth(&self, x: &f32, y: &f32) -> f32;
}

/// Read data from test_bathy_3.nc netcdf3 file that contains x, y, and depth
mod cartesian {

    use std::path::Path;
    use netcdf3::FileReader;
    use crate::{interpolator, bathymetry::BathymetryData};

    /// A struct that stores a netcdf3 dataset named test_bathy_3.nc with methods to access and find nearest values
    pub(crate) struct CartesianFile {
        variables: (Vec<f32>, Vec<f32>, Vec<f32>),
    }

    impl BathymetryData for CartesianFile {
        fn get_depth(&self, x: &f32, y: &f32) -> f32 {
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
        /// `usize` : index of closest value
        /// 
        /// # Panics
        /// This function will panic if direction is not either `"x"` or `"y"`
        fn nearest(&self, target: &f32, direction: &str) -> usize {
            let arr = match direction {
                "x" => &self.variables.0,
                "y" => &self.variables.1,
                _ => todo!("Input a valid option"),
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
        /// `(usize, usize)`: a tuple of x index and y index
        /// 
        /// # Panics
        /// This function will never panic, but if given an out of bounds point,
        /// it will return the closest edge.
        fn nearest_point(&self, x: &f32, y: &f32) -> (usize, usize) {
            let indy = self.nearest(y, "y");
            let indx = self.nearest(x, "x");
            (indx, indy)
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
        /// `Vec<(usize, usize)>` : indices for the four corners surrounding the
        /// given indices in clockwise order.
        /// 
        /// # Panics
        /// This function will not panic, but be aware that it can return values
        /// that are out of bounds to the array.
        fn four_corners(&self, indx: &usize, indy: &usize) -> Vec<(usize, usize)> {
            // clockwise in order
            vec![
                (*indx, indy+1),
                (indx+1, *indy),
                (*indx, indy-1),
                (indx-1, *indy)
            ]
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
        /// `f32` : the depth at the target point
        /// 
        /// # Panics
        /// This function will panic if a point in points is out of bounds of the
        /// array. This function will also panic if `interpolator::bilinear` panics.
        fn interpolate(&self, points: &Vec<(usize, usize)>, target: &(f32, f32)) -> f32 {
            let pts = vec![
                (points[0].0 as i32, points[0].1 as i32, self.depth_from_arr(&points[0].0, &points[0].1)),
                (points[1].0 as i32, points[1].1 as i32, self.depth_from_arr(&points[1].0, &points[1].1)),
                (points[2].0 as i32, points[2].1 as i32, self.depth_from_arr(&points[2].0, &points[2].1)),
                (points[3].0 as i32, points[3].1 as i32, self.depth_from_arr(&points[3].0, &points[3].1)),
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
        /// `f32` : depth
        /// 
        /// # Panics
        /// This function will panic if combined index is out of bounds.
        fn depth_from_arr(&self, indx: &usize, indy: &usize) -> f32 {
            let index = self.variables.0.len() * indy + indx;
            self.variables.2[index]
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
        /// `f32` : depth at the point
        /// 
        /// # Panics
        /// This will panic if the functions `nearest_point`, `four_corners`, or
        /// `interpolate` panic.
        fn depth(&self, x: &f32, y: &f32) -> f32 {
            let nearest_pt = self.nearest_point(x, y);
            let edge_points = self.four_corners(&nearest_pt.0, &nearest_pt.1);
            let depth = self.interpolate(&edge_points, &(*x, *y));
            depth
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
        dbg!(data.nearest(&5499.0, "x"));
        assert!(data.nearest(&5499.0, "x") == 11);
    }

    #[test]
    // view the output from four_corners function
    fn test_get_corners() {
        let data = CartesianFile::new(Path::new("data/test_bathy_3.nc"));
        dbg!(data.four_corners(&10, &10));
    }

    #[test]
    // check values inside the four quadrants but not on grid points
    fn test_depth() {
        let data = CartesianFile::new(Path::new("data/test_bathy_3.nc"));
        assert!(data.get_depth(&10099.0, &5099.0) == 20.0);
        assert!(data.get_depth(&30099.0, &5090.0) == 5.0);
        assert!(data.get_depth(&10099.0, &15099.0) == 10.0);
        assert!(data.get_depth(&30099.0, &15099.0) == 15.0);
    }

}