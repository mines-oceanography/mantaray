//! Structure for containing the results of ray tracing a bundle of rays.

use ode_solvers::dop_shared::SolverResult;
use serde::{Deserialize, Serialize};

use crate::{
    ray_result::RayResult,
    wave_ray_path::{State, Time},
    write_json::WriteJson,
};

#[derive(Serialize, Deserialize, PartialEq, Debug)]
/// Structure for containing the results of ray tracing a bundle of rays.
/// Derives `serde` traits for serialization and deserialization.
///
/// # Note:
/// There is no constructor for this struct, as it is created from a
/// `Vec<Option<SolverResult<Time, State>>>` using the `From` trait.
pub struct BundleResult {
    /// a vector containing a `RayResults` for each ray in the bundle
    rays: Vec<RayResult>,
}

impl WriteJson for BundleResult {}

impl From<Vec<Option<SolverResult<Time, State>>>> for BundleResult {
    /// convert from `Vec<Option<SolverResult<Time, State>>>` (the output from trace_many) to `BundleResult`
    fn from(value: Vec<Option<SolverResult<Time, State>>>) -> Self {
        let mut rays = Vec::new();

        for solver_result in value.iter().flatten() {
            rays.push(solver_result.clone().into());
        }

        BundleResult { rays }
    }
}

#[cfg(test)]
mod test_ray_bundle {
    use super::*;

    #[test]
    fn test_as_json() {
        let ray1: Option<SolverResult<Time, State>> = Some(SolverResult::default());
        let ray2 = Some(SolverResult::default());
        let ray3 = None;

        let result = vec![ray1, ray2, ray3];

        let bundle_result = BundleResult::from(result);

        assert_eq!(bundle_result.as_json(), "{\"rays\":[{\"t_vec\":[],\"x_vec\":[],\"y_vec\":[],\"kx_vec\":[],\"ky_vec\":[]},{\"t_vec\":[],\"x_vec\":[],\"y_vec\":[],\"kx_vec\":[],\"ky_vec\":[]}]}");
    }
}
