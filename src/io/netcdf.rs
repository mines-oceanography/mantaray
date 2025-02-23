use std::collections::HashMap;

use tracing::trace;

use super::{Dataset, LinearFit};
use crate::bathymetry::BathymetryData;
use crate::datatype::{Gradient, Point};
use crate::error::{Error, Result};

/// Implement the Dataset trait for the netcdf::File
impl Dataset for netcdf::File {
    fn dimension_len(&self, name: &str) -> Result<usize> {
        Ok(self.dimension_len(name).unwrap())
    }

    fn values(&self, name: &str) -> Result<ndarray::ArrayD<f64>> {
        Ok(self.variable(name).unwrap().get::<f64, _>(..).unwrap())
    }

    // Missing get full variable (such as all x values), and get size.
    fn get_variable(&self, name: &str, i: usize, j: usize) -> Result<f32> {
        Ok(self
            .variable(name)
            .unwrap()
            .get_value::<f32, _>([i, j])
            .unwrap())
    }
}

struct RegularGrid {
    dataset: netcdf::File,
    x_size: usize,
    x_map: LinearFit<f64>,
    y_size: usize,
    y_map: LinearFit<f64>,
    // Save the dimensions order: ij or ji
    dimension_order: HashMap<String, String>,
}

impl RegularGrid {
    /*
        fn validate_as_regular_grid() {
            // Open a full 1D array for x and another for y
            // calculate delta and confirm that all are <0.01 diference (criteria for linear)
            //
        }
    */
    #[allow(dead_code)]
    fn open<P>(file: P, varname_x: &str, varname_y: &str) -> Self
    where
        P: AsRef<std::path::Path>,
    {
        // confirm it is linear
        // Define A & B coefficients
        // get i_size and j_size

        let dataset = netcdf::open(&file).unwrap();

        // Identify the variables that have the user defined dimensions
        // and create a map on the dimenson order
        let varnames = &dataset
            .variables()
            .into_iter()
            .filter(|v| {
                v.dimensions()
                    .into_iter()
                    .map(|v| v.name() == varname_x)
                    .any(|v| v)
            })
            .filter(|v| {
                v.dimensions()
                    .into_iter()
                    .map(|v| v.name() == varname_y)
                    .any(|v| v)
            })
            .filter_map(|v| {
                match &v
                    .dimensions()
                    .into_iter()
                    .map(|v| v.name())
                    .collect::<Vec<_>>()[..]
                {
                    [varname_x, varname_y] => Some((v.name(), "xy".to_string())),
                    [varname_y, varname_x] => Some((v.name(), "yx".to_string())),
                    _ => None,
                }
            })
            .collect::<Vec<_>>();

        let dimension_order: HashMap<String, String> = HashMap::from_iter(varnames.iter().cloned());

        let x_size = dataset.dimension_len(varname_x).unwrap();

        let x_map = LinearFit::from_fit(dataset.values(varname_x).unwrap()).unwrap();

        let y_size = dataset.dimension_len(varname_y).unwrap();
        let y_map = LinearFit::from_fit(dataset.values(varname_x).unwrap()).unwrap();

        Self {
            dataset,
            x_size,
            x_map,
            y_size,
            y_map,
            dimension_order,
        }
    }

    #[allow(dead_code)]
    /// Get the nearest `varname` value to the given `x` and `y` coordinates
    fn nearest(&self, varname: &str, point: Point<f64>) -> Result<f32> {
        let i = self.x_map.predict(*point.x()).round() as usize;
        if i >= self.x_size {
            return Err(Error::IndexOutOfBounds);
        }
        let j = self.y_map.predict(*point.y()).round() as usize;
        if j >= self.y_size {
            return Err(Error::IndexOutOfBounds);
        }

        match self.dimension_order.get(varname) {
            Some(v) => match v.as_str() {
                "xy" => {
                    trace!("Assuming dimension order is 'xy'");
                    self.dataset.get_variable(varname, i, j)
                }
                "yx" => {
                    trace!("Assuming dimension order is 'yx'");
                    self.dataset.get_variable(varname, j, i)
                }
                _ => return Err(Error::Undefined),
            },
            _ => return Err(Error::Undefined),
        }
    }
}

#[allow(dead_code)]
impl BathymetryFromNetCDF {
    pub fn new<P>(file: P, x_name: &str, y_name: &str, depth_name: String) -> Self
    where
        P: AsRef<std::path::Path>,
    {
        let file = netcdf::open(&file).unwrap();

        let x: Vec<f32> = file
            .variable(x_name)
            .expect("Could not find variable 'x'")
            .get::<f32, _>(..)
            .expect("Could not get value of variable 'x'")
            .into_raw_vec();

        let y: Vec<f32> = file
            .variable(y_name)
            .expect("Could not find variable 'y'")
            .get::<f32, _>(..)
            .expect("Could not get value of variable 'y'")
            .into_raw_vec();

        Self {
            file,
            x,
            y,
            depth_name,
        }
    }
}

impl BathymetryFromNetCDF {
    fn nearest_location_index(&self, x0: &f32, y0: &f32) -> Result<(usize, usize)> {
        let i = self
            .x
            .iter()
            .enumerate()
            .map(|v| (v.0, (x0 - *v.1).abs()))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .expect("Could not find minimum")
            .0;

        let j = self
            .y
            .iter()
            .enumerate()
            .map(|v| (v.0, (y0 - *v.1 as f32).abs()))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .expect("Could not find minimum")
            .0;

        Ok((i, j))
    }

    fn depth_by_index(&self, i: usize, j: usize) -> f32 {
        let z0 = self
            .file
            .variable(&self.depth_name)
            .expect("Could not find variable 'depth'")
            .get_value::<f32, _>([j, i])
            .expect("Could not get value of variable 'depth'");

        z0
    }
}

impl BathymetryData for BathymetryFromNetCDF {
    fn depth(&self, point: &Point<f32>) -> Result<f32> {
        let (i, j) = self.nearest_location_index(point.x(), point.y())?;
        Ok(self.depth_by_index(i, j))
    }

    fn depth_and_gradient(&self, point: &Point<f32>) -> Result<(f32, Gradient<f32>)> {
        let (i, j) = self.nearest_location_index(point.x(), point.y())?;
        let z0 = self.depth_by_index(i, j);

        let delta_2x = self.x[i + 1] - self.x[i - 1];
        let dzdx = (self.depth_by_index(i + 1, j) - self.depth_by_index(i - 1, j)) / delta_2x;

        let delta_2y = self.y[j + 1] - self.y[j - 1];
        let dzdy = (self.depth_by_index(i, j + 1) + self.depth_by_index(i, j - 1)) / delta_2y;

        Ok((z0, Gradient::new(dzdx, dzdy)))
    }
}

/*
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_depth() {
        let bathymetry =
            BathymetryFromNetCDF::new("gaussian_island.nc", &"x", &"y", "depth".to_string());

        let depth = bathymetry.depth(&800.0, &128.0).unwrap();
        assert_eq!(depth, 21.388426);
        bathymetry
            .depth_and_gradient(&1_000.0, &3_141.0)
            .unwrap();
    }

    #[test]
    fn read_depth_and_gradient() {
        let bathymetry =
            BathymetryFromNetCDF::new("gaussian_island.nc", &"x", &"y", "depth".to_string());

        let (depth, (dx, dy)) = bathymetry.depth_and_gradient(&800.0, &128.0).unwrap();
        assert_eq!(depth, 21.388426);
        assert_eq!(dx, 0.056151398);
        assert_eq!(dy, 0.09486287);
    }
}
*/
