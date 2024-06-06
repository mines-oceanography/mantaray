//! RayResults struct which holds the results of the ray tracing as vectors.
//! Contains methods to convert from `SolverResult` to `RayResults`.

use ode_solvers::dop_shared::SolverResult;
use serde::Deserialize;
use serde::Serialize;

use crate::wave_ray_path::{State, Time};
use crate::write_json::WriteJson;

#[derive(Serialize, Deserialize, PartialEq, Debug)]
/// struct to hold the results of the ray tracing simulation as vectors. Note
/// that the vectors are not indexed by time, but by the number of steps of the
/// simulation.
pub struct RayResult {
    /// vector of time values
    t: Vec<f64>,
    /// vector of x location values
    x: Vec<f64>,
    /// vector of y location values
    y: Vec<f64>,
    /// vector of kx values
    kx: Vec<f64>,
    /// vector of ky values
    ky: Vec<f64>,
}

#[allow(dead_code)]
impl RayResult {
    /// Create a new RayResults struct with the given vectors.
    ///
    /// # Arguments
    ///
    /// `t` : `Vec<f64>`
    /// - a vector of time values
    ///
    /// `x` : `Vec<f64>`
    /// - a vector of x values
    ///
    /// `y` : `Vec<f64>`
    /// - a vector of y values
    ///
    /// `kx` : `Vec<f64>`
    /// - a vector of kx values
    ///
    /// `ky` : `Vec<f64>`
    /// - a vector of ky values
    ///
    /// # Returns
    ///
    /// constructed `RayResults` struct
    pub fn new(t: Vec<f64>, x: Vec<f64>, y: Vec<f64>, kx: Vec<f64>, ky: Vec<f64>) -> Self {
        RayResult { t, x, y, kx, ky }
    }
}

impl WriteJson for RayResult {}

impl From<SolverResult<Time, State>> for RayResult {
    /// convert the SolverResult to a RayResults struct
    fn from(value: SolverResult<Time, State>) -> Self {
        let (x_out, y_out) = value.get();

        let mut t: Vec<f64> = vec![];
        let mut x: Vec<f64> = vec![];
        let mut y: Vec<f64> = vec![];
        let mut kx: Vec<f64> = vec![];
        let mut ky: Vec<f64> = vec![];

        for (i, _) in x_out.iter().enumerate() {
            if y_out[i][0].is_nan()
                || y_out[i][1].is_nan()
                || y_out[i][2].is_nan()
                || y_out[i][3].is_nan()
            {
                break;
            }
            t.push(x_out[i]);
            x.push(y_out[i][0]);
            y.push(y_out[i][1]);
            kx.push(y_out[i][2]);
            ky.push(y_out[i][3]);
        }

        RayResult::new(t, x, y, kx, ky)
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
    /// test the to_json_string method
    fn test_to_json_string() {
        let ray_results = RayResult::new(vec![1.0], vec![2.0], vec![3.0], vec![4.0], vec![5.0]);

        let json_string = ray_results.to_json_string();

        assert_eq!(
            json_string,
            "{\"t\":[1.0],\"x\":[2.0],\"y\":[3.0],\"kx\":[4.0],\"ky\":[5.0]}"
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

        let json_string = rr.to_json_string();

        assert_eq!(
            json_string,
            "{\"t\":[0.0],\"x\":[1.0],\"y\":[1.0],\"kx\":[1.0],\"ky\":[1.0]}"
        );
    }
}
