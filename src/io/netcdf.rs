use crate::bathymetry::BathymetryData;
use crate::datatype::{Gradient, Point};
use crate::error::{Error, Result};

#[allow(dead_code)]
#[derive(Debug)]
/// Holds and apply a linear relationship
///
/// This was originally created to handle the dimensions of regular girds,
/// such latitude and longitude in a regularly spaced dataset, providing a
/// cheap conversion between the dimension, such as latitute, to the
/// correspondent index position.
struct LinearFit<T> {
    slope: T,
    intercept: T,
}

impl<T> LinearFit<T>
where
    T: Copy,
    T: std::ops::Sub<Output = T>,
    T: std::ops::Mul<Output = T>,
{
    #[allow(dead_code)]
    /// Predict the index position of a given value
    fn predict(&self, x: T) -> T {
        (x - self.intercept) * self.slope
    }
}

impl LinearFit<f64> {
    /// Create a new LinearFit from a vector of values
    ///
    /// For a given vector of values, calculates the best linear relationship
    /// between the values and its index position with the purpose to estimate
    /// the index position of the closest value to a given target.
    ///
    /// This procedure also validates if a linear relationship is a good
    /// approximation with a threshold of 0.5% of tolerance.
    fn from_fit(x: ndarray::ArrayD<f64>) -> Result<LinearFit<f64>> {
        let dx = &x.slice(ndarray::s![1..]) - &x.slice(ndarray::s![..-1]);
        let slope = dx.mean().unwrap();
        let criteria = ((dx - slope) / slope)
            .abs()
            .into_iter()
            .map(|v| v > 0.005)
            .any(|v| v);
        if criteria {
            return Err(Error::Undefined);
        }
        Ok(LinearFit {
            slope,
            intercept: x[0],
        })
    }

    /*
    fn from_fit(x: Vec<f32>) -> Result<Self> {
        let delta: Vec<_> = x.windows(2).map(|v| v[1] - v[0]).collect();
        let slope = (delta.len() as f32 - 1.0) / delta.iter().sum::<f32>();
        let threshold = delta
            .iter()
            .map(|v| (v - slope).abs() / slope)
            .any(|v| v > 0.005);
        if threshold {
            return Err(Error::Undefined);
        }

        let intercept = x[0];
        Ok(LinearFit { slope, intercept })
    }
    */
}

trait Dataset {
    fn dimension_len(&self, name: &str) -> Result<usize>;
    fn values(&self, name: &str) -> Result<ndarray::ArrayD<f64>>;
    fn get_variable(&self, name: &str, i: usize, j: usize) -> Result<f32>;
}

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
    x: LinearFit<f32>,
    x_size: usize,
    y: LinearFit<f32>,
    y_size: usize,
}

impl RegularGrid {
    /*
        fn validate_as_regular_grid() {
            // Open a full 1D array for x and another for y
            // calculate delta and confirm that all are <0.01 diference (criteria for linear)
            //
        }
    */
    fn open<P>(file: P) -> Self
    where
        P: AsRef<std::path::Path>,
    {
        // confirm it is linear
        // Define A & B coefficients
        // get x_size and y_size

        let dataset = netcdf::open(&file).unwrap();
        let varname_x = "ETOPO05_X";
        let varname_y = "ETOPO05_Y";

        let x_size = dataset.dimension(varname_x).unwrap().len();
        let x = dataset
            .variable(varname_x)
            .unwrap()
            .get::<f64, _>(..)
            .unwrap();
        let map_i = LinearFit::from_fit(x.iter().map(|v| *v as f32).collect::<Vec<_>>()).unwrap();
        let y_size = dataset.dimension(varname_y).unwrap().len();
        let y = dataset
            .variable(varname_y)
            .unwrap()
            .get::<f64, _>(..)
            .unwrap();
        let map_j = LinearFit::from_fit(y.iter().map(|v| *v as f32).collect::<Vec<_>>()).unwrap();

        Self {
            dataset,
            x: map_i,
            x_size,
            y: map_j,
            y_size,
        }
    }

    /// Get the nearest `varname` value to the given `x` and `y` coordinates
    fn nearest(&self, varname: &str, x: f32, y: f32) -> Result<f32> {
        let j = self.x.predict(x).round() as usize;
        let i = self.y.predict(y).round() as usize;
        let z = self.dataset.get_variable(varname, i, j).unwrap();
        Ok(z)
    }
}

pub(crate) struct BathymetryFromNetCDF {
    file: netcdf::File,
    x: Vec<f32>,
    y: Vec<f32>,
    depth_name: String,
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
