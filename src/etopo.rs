//! Module to read a netcdf file with bathymetry data.
//! 
//! The module is currently only tested for test_bathy_3.nc.
//! 
//! Requires netcdf3 crate. Will be using a interpolation crate in the future.

mod etopo {

    use std::path::Path;
    use netcdf3::{FileReader, DataSet, DataVector, DataType, Version, DimensionType};

    /// A struct that stores a netcdf dataset with methods to access and find nearest values
    pub(crate) struct BathyData {
        path: String,
        xname: String,
        yname: String,
        depth_name: String,
        data: FileReader,
        variables: (Vec<f32>, Vec<f32>, Vec<f64>),
    }

    impl BathyData {
        /// Construct BathyData
        /// 
        /// # Arguments
        /// `path` : `String`
        /// - a path to the location of the netcdf file
        /// 
        /// `xname` : `String`
        /// - name of the x variable in the netcdf file
        /// 
        /// `yname` : `String`
        /// - name of the y variable in the netcdf file
        /// 
        /// `depth_name` : `String`
        /// - name of the depth variable in the netcdf file
        /// 
        /// # Returns
        /// `Self` : an initialized BathyData
        /// 
        /// # Panics
        /// `new` will panic if the data type is invalid or if any of the names are invalid.
        fn new(path: String, xname: String, yname: String, depth_name: String) -> Self {
            let mut data = FileReader::open(Path::new(&path)).unwrap();
            let variables = (
                data.read_var_f32(&xname).unwrap(),
                data.read_var_f32(&yname).unwrap(),
                data.read_var_f64(&depth_name).unwrap()
            );
            BathyData { path, xname, yname, depth_name, data, variables }
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

            let target = target as f32;

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
        /// Get four adjecent points
        /// 
        /// 
        fn four_corners() -> [(f64, f64); 4] {
            todo!()
        }
        /// Interpolate the depth
        /// 
        /// 
        fn interpolate() -> f64 {
            todo!()
        }
        /// Return the depth at x, y
        fn depth(x: f64, y: f64) -> f64 {
            todo!()
        }

    }

    /// this function creates a pointer to the struct and returns it.
    pub(crate) fn test_bathy_3_data() -> Box<BathyData> {
        let path = String::from("C:\\Users\\bairv\\ray_tracing_etopo5\\data\\test_bathy_3.nc");
        let xname = String::from("x");
        let yname =String::from("y");
        let depth_name = String::from("depth");
        Box::new(BathyData::new(path, xname, yname, depth_name))
    }

    /// a function to open the etopo5.nc file and return pointers to variables
    pub(crate) fn open_variables() -> (Vec<f32>, Vec<f32>, Vec<f64>) {
        let depth_data = test_bathy_3_data();

        depth_data.variables
    }

    /// a function to view the file data
    pub(crate) fn view_file_data() {
        let mut file_reader: FileReader = FileReader::open(r"C:\Users\bairv\ray_tracing_etopo5\data\test_bathy_3.nc").unwrap();

        for var in file_reader.data_set().get_vars() {
            println!("{}", var.name());
        }
    }

    /// this function creates the dataset and calls the nearest function
    pub(crate) fn get_nearest() -> usize {
        let depth_data = test_bathy_3_data();

        depth_data.nearest(5499.0, "x")
    }
}

#[cfg(test)]
mod test_netcdf {

    use crate::etopo::etopo::get_nearest;
    use super::etopo::{open_variables, view_file_data};

    #[test]
    /// view the variable names in the file
    fn test_view_file_data() {
        view_file_data();
    }

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

}