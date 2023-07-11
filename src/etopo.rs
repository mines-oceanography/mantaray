//! Module to read a netcdf file with bathymetry data.
//! 
//! The module is currently only tested for etopo5.nc from oceansdb python library.
//! 

mod etopo {

    use std::{collections::HashMap,};
    use netcdf3::{FileReader, DataSet, DataVector, DataType, Version, DimensionType};

    /// find the closest value to a given latitude or longitude using a simple
    /// linear search. This is very slow and takes linear time, where using a tree approach could make it run much faster.
    pub(crate) fn closest_value_index(target: f64, arr: &Vec<f64>) -> usize {
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

    /// a function to open the etopo5.nc file and return pointers to variables
    pub(crate) fn open_variables() -> (Box<Vec<f64>>, Box<Vec<f64>>, Box<Vec<f32>>) {
        let mut file_reader: FileReader = FileReader::open(r"C:\Users\bairv\ray_tracing_etopo5\src\etopo5.nc").unwrap();

        let etopo05_y = file_reader.read_var_f64("ETOPO05_Y").unwrap();
        let etopo05_x = file_reader.read_var_f64("ETOPO05_X").unwrap();
        let rose = file_reader.read_var_f32("ROSE").unwrap();

        (Box::new(etopo05_y), Box::new(etopo05_x), Box::new(rose))
    }

}

#[cfg(test)]
mod test_netcdf {

    use std::{collections::HashMap};
    use netcdf3::{FileReader, DataSet, DataVector, DataType, Version, DimensionType};

    use super::etopo::{closest_value_index, open_variables};

    #[test]
    /// test access to variables created by open_variables
    fn test_open_variables() {
        let (lat, lon, depth) = open_variables();
        println!("{}, {}, {}", lat[0], lon[0], depth[0]);
    }

    #[test]
    /// This function will get the latitude, longitude, and depth from the file.
    /// Then, it will use the closest_value_index to find the closest values.
    /// These should be the exact same because the inputs are taken from the
    /// file.
    fn test_closest_index() {
        let (etopo05_y, etopo05_x, rose) = open_variables();

        let i = 2000;
        let j = 100;
        let target_lat = etopo05_y[i];
        let target_lon = etopo05_x[j];
        let test_index = i * etopo05_x.len() + j;


        let res_lat_i = closest_value_index(target_lat, &etopo05_y);
        let res_lon_i = closest_value_index(target_lon, &etopo05_x);

        let depth_index = res_lat_i * etopo05_x.len() + res_lon_i;
        println!("{depth_index}");

        println!("latitude: {}, longitude: {}, depth: {}", etopo05_y[res_lat_i], etopo05_x[res_lon_i], rose[depth_index]);

        assert_eq!(depth_index, test_index)
    }
}