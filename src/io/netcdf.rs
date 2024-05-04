use crate::bathymetry::BathymetryData;
use crate::error::Error;

pub(crate) struct BathymetryFromNetCDF {
    file: netcdf::File,
    x_name: String,
    y_name: String,
    depth_name: String,
}

impl BathymetryFromNetCDF {
    pub fn new<P>(file: P, x_name: String, y_name: String, depth_name: String) -> Self
    where
        P: AsRef<std::path::Path>,
    {
        Self {
            file: netcdf::open(&file).unwrap(),
            x_name,
            y_name,
            depth_name,
        }
    }
}

impl BathymetryData for BathymetryFromNetCDF {
    fn get_depth(&self, x0: &f32, y0: &f32) -> Result<f32, Error> {
        // Load x
        let x = self
            .file
            .variable(&self.x_name)
            .expect("Could not find variable 'x'")
            .get::<f64, _>(..)
            .expect("Could not get value of variable 'x'");

        // Find the closest value to x0
        let i = x
            .iter()
            .enumerate()
            .map(|v| (v.0, (x0 - *v.1 as f32).abs()))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .expect("Could not find minimum")
            .0;

        // Load y
        let y = self
            .file
            .variable(&self.y_name)
            .expect("Could not find variable 'y'")
            .get::<f64, _>(..)
            .expect("Could not get value of variable 'y'");

        // Find the closest value to y0
        let j = y
            .iter()
            .enumerate()
            .map(|v| (v.0, (y0 - *v.1 as f32).abs()))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .expect("Could not find minimum")
            .0;

        // Load z at (i, j)
        let depth = self
            .file
            .variable(&self.depth_name)
            .expect("Could not find variable 'depth'")
            .get_value::<f64, _>([j, i])
            .expect("Could not get value of variable 'depth'");

        Ok(depth as f32)
    }

    fn get_depth_and_gradient(&self, x: &f32, y: &f32) -> Result<(f32, (f32, f32)), Error> {
        Ok((0.0, (0.0, 0.0)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create() {
        let bathymetry = BathymetryFromNetCDF::new(
            "bathymetry.nc",
            "ETOPO05_X".to_string(),
            "ETOPO05_Y".to_string(),
            "ROSE".to_string(),
        );

        let z = bathymetry.get_depth(&11.234, &0.0).unwrap();
        dbg!(z);
        assert!(false);
    }

    /*
    #[test]
    fn test_bathymetry_from_netcdf() {
        let bathymetry = BathymetryFromNetCDF::new(
            "data/bathymetry.nc",
            "x".to_string(),
            "y".to_string(),
            "depth".to_string(),
        );
        assert_eq!(bathymetry.get_depth(&0.0, &0.0), &0.0);
        assert_eq!(
            bathymetry.get_depth_and_gradient(&0.0, &0.0),
            (0.0, (0.0, 0.0))
        );
    }
    */
}
