//! RayResults struct which holds the results of the ray tracing as vectors.
//! Contains methods to convert from `SolverResult` and to `RayResults` and
//! implements the `Writable` trait. 

use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

use ode_solvers::dop_shared::SolverResult;

use crate::error::Error;
use crate::wave_ray_path::{State, Time};
use crate::writable::Writable;

#[allow(dead_code)]
/// struct to hold the results of the ray tracing simulation as vectors. Note
/// that the vectors are not indexed by time, but by the number of steps of the
/// simulation.
pub struct RayResults {
    t_vec: Vec<f64>,
    x_vec: Vec<f64>,
    y_vec: Vec<f64>,
    kx_vec: Vec<f64>,
    ky_vec: Vec<f64>,
}

impl RayResults {
    /// Create a new RayResults struct with the given vectors.
    pub fn new(
        t_vec: Vec<f64>,
        x_vec: Vec<f64>,
        y_vec: Vec<f64>,
        kx_vec: Vec<f64>,
        ky_vec: Vec<f64>,
    ) -> Self {
        RayResults {
            t_vec,
            x_vec,
            y_vec,
            kx_vec,
            ky_vec,
        }
    }
}

impl From<SolverResult<Time, State>> for RayResults {
    /// convert the SolverResult to a RayResults struct
    fn from(value: SolverResult<Time, State>) -> Self {
        let (x_out, y_out) = value.get();

        let mut t_vector = vec![];
        let mut x_vector: Vec<f64> = vec![];
        let mut y_vector: Vec<f64> = vec![];
        let mut kx_vector: Vec<f64> = vec![];
        let mut ky_vector: Vec<f64> = vec![];

        for (i, _) in x_out.iter().enumerate() {
            if y_out[i][0].is_nan() {
                break;
            }
            t_vector.push(x_out[i]);
            x_vector.push(y_out[i][0]);
            y_vector.push(y_out[i][1]);
            kx_vector.push(y_out[i][2]);
            ky_vector.push(y_out[i][3]);
        }

        RayResults::new(t_vector, x_vector, y_vector, kx_vector, ky_vector)
    }
}

impl Writable for RayResults {
    fn write_to_csv_file(&self, path: &std::path::Path, delimiter: &str) -> Result<(), Error> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        writeln!(
            &mut writer,
            "t{}x{}y{}kx{}ky{}",
            delimiter, delimiter, delimiter, delimiter, delimiter
        )?;
        for i in 0..self.t_vec.len() {
            writeln!(
                &mut writer,
                "{}{}{}{}{}{}{}{}{}{}",
                self.t_vec[i],
                delimiter,
                self.x_vec[i],
                delimiter,
                self.y_vec[i],
                delimiter,
                self.kx_vec[i],
                delimiter,
                self.ky_vec[i],
                delimiter
            )?;
        }

        writer.flush()?;
        Ok(())
    }
}
