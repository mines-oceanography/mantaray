//! RayResults struct which holds the results of the ray tracing as vectors.
//! Contains methods to convert from `SolverResult` and to `RayResults` and
//! write using serde and serde_json.

use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

use ode_solvers::dop_shared::SolverResult;
use serde::Deserialize;
use serde::Serialize;

use crate::error::Result;
use crate::wave_ray_path::{State, Time};

#[derive(Serialize, Deserialize, PartialEq, Debug)]
/// struct to hold the results of the ray tracing simulation as vectors. Note
/// that the vectors are not indexed by time, but by the number of steps of the
/// simulation.
pub(crate) struct RayResult {
    t_vec: Vec<f64>,
    x_vec: Vec<f64>,
    y_vec: Vec<f64>,
    kx_vec: Vec<f64>,
    ky_vec: Vec<f64>,
}

#[allow(dead_code)]
impl RayResult {
    /// Create a new RayResults struct with the given vectors.
    ///
    /// # Arguments
    ///
    /// `t_vec` : `Vec<f64>`
    /// - a vector of time values
    ///
    /// `x_vec` : `Vec<f64>`
    /// - a vector of x values
    ///
    /// `y_vec` : `Vec<f64>`
    /// - a vector of y values
    ///
    /// `kx_vec` : `Vec<f64>`
    /// - a vector of kx values
    ///
    /// `ky_vec` : `Vec<f64>`
    /// - a vector of ky values
    ///
    /// # Returns
    ///
    /// constructed `RayResults` struct
    pub(crate) fn new(
        t_vec: Vec<f64>,
        x_vec: Vec<f64>,
        y_vec: Vec<f64>,
        kx_vec: Vec<f64>,
        ky_vec: Vec<f64>,
    ) -> Self {
        RayResult {
            t_vec,
            x_vec,
            y_vec,
            kx_vec,
            ky_vec,
        }
    }

    /// Convert the `RayResults` struct to a JSON string.
    ///
    /// # Returns
    ///
    /// JSON string of the `RayResults` struct
    pub(crate) fn as_json(&self) -> String {
        serde_json::to_string(&self).unwrap()
    }

    /// Write the `RayResults` struct to a writer.
    ///
    /// # Arguments
    ///
    /// `writer` : `&mut W`
    /// - object that implements `Write` to write the `RayResults` struct to
    ///
    /// # Returns
    ///
    /// `Ok(usize)` : the number of bytes written
    ///
    /// `Err(Error)` : an error occurred while writing
    ///
    /// # Note
    ///
    /// This method writes the `RayResults` struct as a JSON string.
    pub(crate) fn write<W: Write>(&self, writer: &mut W) -> Result<usize> {
        writer.write_all(self.as_json().as_bytes())?;
        writer.flush()?;
        Ok(self.as_json().as_bytes().len())
    }

    /// Save the `RayResults` struct to a file at the given path.
    ///
    /// # Arguments
    ///
    /// `path` : `&Path`
    /// - the path to save the `RayResults` struct to
    ///
    /// # Returns
    ///
    /// `Ok(usize)` : the number of bytes written
    ///
    /// `Err(Error)` : an error occurred while writing
    ///
    /// # Note
    ///
    /// This method writes the `RayResults` struct as a JSON string at the given file path.
    pub(crate) fn save_file(&self, path: &Path) -> Result<usize> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        self.write(&mut writer)
    }
}

impl From<SolverResult<Time, State>> for RayResult {
    /// convert the SolverResult to a RayResults struct
    fn from(value: SolverResult<Time, State>) -> Self {
        let (x_out, y_out) = value.get();

        let mut t_vector = vec![];
        let mut x_vector: Vec<f64> = vec![];
        let mut y_vector: Vec<f64> = vec![];
        let mut kx_vector: Vec<f64> = vec![];
        let mut ky_vector: Vec<f64> = vec![];

        for (i, _) in x_out.iter().enumerate() {
            if y_out[i][0].is_nan()
                || y_out[i][1].is_nan()
                || y_out[i][2].is_nan()
                || y_out[i][3].is_nan()
            {
                break;
            }
            t_vector.push(x_out[i]);
            x_vector.push(y_out[i][0]);
            y_vector.push(y_out[i][1]);
            kx_vector.push(y_out[i][2]);
            ky_vector.push(y_out[i][3]);
        }

        RayResult::new(t_vector, x_vector, y_vector, kx_vector, ky_vector)
    }
}

#[cfg(test)]
mod test_ray_result {

    use super::*;

    #[test]
    /// test the converted RayResults struct from a SolverResult with constructor
    fn test_ray_result() {
        let solver_result: SolverResult<Time, State> = SolverResult::default();

        let converted_ray_results = RayResult::from(solver_result);

        let constructed_ray_results = RayResult::new(vec![], vec![], vec![], vec![], vec![]);

        assert_eq!(converted_ray_results, constructed_ray_results);
    }

    #[test]
    /// test the as_json method
    fn test_as_json() {
        let ray_results = RayResult::new(vec![1.0], vec![2.0], vec![3.0], vec![4.0], vec![5.0]);

        let json_string = ray_results.as_json();

        assert_eq!(
            json_string,
            "{\"t_vec\":[1.0],\"x_vec\":[2.0],\"y_vec\":[3.0],\"kx_vec\":[4.0],\"ky_vec\":[5.0]}"
        );
    }

    #[test]
    /// test NaN. when converting the `SolverResult` to `RayResult`, if an entry
    /// in the `SolverResult` has a NaN value, then that value and all after it
    /// are not included in the converted `RayResult`.
    fn test_nan_ray_result() {
        let sr: SolverResult<Time, State> = SolverResult::new(
            vec![0.0, 1.0, 2.0, 3.0],
            vec![
                State::new(1.0, 1.0, 1.0, 1.0),
                State::new(1.0, f64::NAN, f64::NAN, 1.0),
                State::new(f64::NAN, f64::NAN, f64::NAN, f64::NAN),
                State::new(2.0, 2.0, 2.0, 2.0),
            ],
        );

        let rr: RayResult = sr.into();

        let json_string = rr.as_json();

        assert_eq!(
            json_string,
            "{\"t_vec\":[0.0],\"x_vec\":[1.0],\"y_vec\":[1.0],\"kx_vec\":[1.0],\"ky_vec\":[1.0]}"
        );
    }
}
