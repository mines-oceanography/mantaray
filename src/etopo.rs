//! Module to read a netcdf file with bathymetry data.
//! 
//! The module is currently only tested for etopo5.nc.
//! 
//! Requires netcdf3 crate. Will be using a interpolation crate in the future.

trait BathymetryData {
    /// Returns the nearest depth for the given lat, lon coordinate.
    fn get_depth_nearest(&self, lat: &f64, lon: &f64) -> f64;
}

mod etopo {

    use std::path::Path;
    use netcdf3::FileReader;

    use super::BathymetryData;

    /// A struct that stores a netcdf dataset with methods to access and find nearest values
    pub(crate) struct Etopo5 {
        variables: (Vec<f64>, Vec<f64>, Vec<f32>),
    }

    impl BathymetryData for Etopo5 {
        fn get_depth_nearest(&self, lat: &f64, lon: &f64) -> f64 {
            // get the index position for the nearest
            let (x, y) = self.nearest_point(lat, lon);
            // then get actual bathymetry
            self.depth_from_arr(x, y) as f64
        }
    }

    impl Etopo5 {
        /// Construct Etopo5
        /// 
        /// # Arguments
        /// `path` : `&str`
        /// - a path to the location of the netcdf file
        /// 
        /// # Returns
        /// `Self` : an initialized BathyData
        /// 
        /// # Panics
        /// `new` will panic if the data type is invalid or if any of the names
        /// are invalid. But this should never panic for etopo5.nc
        pub(crate) fn new(path: &str) -> Self {
            let mut data = FileReader::open(Path::new(path)).unwrap();
            let variables = (
                data.read_var_f64("ETOPO05_X").unwrap(),
                data.read_var_f64("ETOPO05_Y").unwrap(),
                data.read_var_f32("ROSE").unwrap()
            );
            Etopo5 { variables }
        }
        /// Find nearest point
        /// 
        /// # Arguments
        /// `target` : `f64`
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
        fn nearest(&self, target: f64, direction: &str) -> usize {
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
        /// Returns the nearest x index, y index point to given lat, lon
        /// 
        /// # Arguments
        /// `lat`: `&f64`
        /// - latitude (y) coordinate in range -90.0 to 90.0
        /// 
        /// `lon`: `&f64`
        /// - longitude (x) coordinate in range 0.0 to 359.92
        /// 
        /// # Returns
        /// `(usize, usize)`: a tuple of x index and y index
        /// 
        /// # Panics
        /// This function will never panic, but if given an out of bounds point,
        /// it will return the closest edge.
        fn nearest_point(&self, lat: &f64, lon: &f64) -> (usize, usize) {
            let indy = self.nearest(*lat, "y");
            let indx = self.nearest(*lon, "x");
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
        /// given indices.
        /// 
        /// # Panics
        /// This function will not panic, but be aware that it can return values
        /// that are out of bounds to the array.
        fn four_corners(&self, indx: usize, indy: usize) -> Vec<(usize, usize)> {
            let mut corners = Vec::new();
            corners.push((indy-1, indx));
            corners.push((indy, indx-1));
            corners.push((indy+1, indx));
            corners.push((indy, indx+1));

            corners
        }
        /// Interpolate the depth
        /// 
        /// 
        fn interpolate(&self, points: Vec<(usize, usize)>) -> f64 {
            todo!()
        }
        /// Access values in flattened array as you would a 2d array
        fn depth_from_arr(&self, indx: usize, indy: usize) -> f32 {
            let index = self.variables.0.len() * indy + indx;
            self.variables.2[index]
        }
        /// Return the depth at x, y
        fn depth(x: f64, y: f64) -> f64 {
            todo!()
        }

    }

    /// this function creates a pointer to the struct and returns it.
    pub(crate) fn test_bathy_3_data() -> Box<Etopo5> {
        Box::new(Etopo5::new("data/etopo5.nc"))
    }

    /// a function to open the etopo5.nc file and return pointers to variables
    pub(crate) fn open_variables() -> (Vec<f64>, Vec<f64>, Vec<f32>) {
        let depth_data = test_bathy_3_data();

        depth_data.variables
    }

    /// this function creates the dataset and calls the nearest function
    pub(crate) fn get_nearest() -> usize {
        let depth_data = test_bathy_3_data();

        depth_data.nearest(5499.0, "x")
    }

    /// this function creates the dataset and returns the four corners around a point
    pub(crate) fn get_corners() -> Vec<(usize, usize)> {
        let depth_data = test_bathy_3_data();
        depth_data.four_corners(10, 10)
    }


}

#[cfg(test)]
mod test_netcdf {

    use crate::etopo::etopo::{get_nearest, get_corners};
    use super::{etopo::{open_variables, Etopo5}, BathymetryData};

    #[test]
    /// test access to variables created by open_variables
    fn test_open_variables() {
        let vars = open_variables();
        println!("{}", vars.2[0])
    }

    #[test]
    fn test_nearest() {
        println!("{}", get_nearest())
    }

    #[test]
    fn test_corners() {
        println!("{:?}", get_corners())
    }

    #[test]
    /// This test works for values in the range of etopo5.nc
    fn nearest_bathymetry() {
        // Mines, Golden, CO 39.7510, 254.7774 @ 1747 meters
        // Atacama observatory, -22.9586, 292.2125 @ 4724 meters
        // Titanic, 41.72583043, 310.05917043 @ -3780 meters
        let lat = 41.72583043;
        let lon = 310.05917043;
        // open etopo5 data set
        let etopo_data = Etopo5::new("data/etopo5.nc");
        // use trait function
        dbg!(&etopo_data.get_depth_nearest(&lat, &lon));
        assert!((etopo_data.get_depth_nearest(&lat, &lon) - -3780.0).abs() < f64::EPSILON)
    }


}