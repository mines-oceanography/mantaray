//! Bathymetry

mod array_depth;
pub mod cartesian_netcdf3;
mod constant_depth;
mod constant_slope;

use crate::error::Error;
#[allow(unused_imports)]
pub(super) use array_depth::ArrayDepth;
#[allow(unused_imports)]
pub(super) use cartesian_netcdf3::CartesianNetCDF3;
#[allow(unused_imports)]
pub(super) use constant_depth::ConstantDepth;
#[allow(unused_imports)]
pub(super) use constant_slope::ConstantSlope;

/// A trait used to give the function get_depth
pub trait BathymetryData: Sync {
    /// Returns the nearest depth for the given x, y coordinate.
    fn depth(&self, x: &f32, y: &f32) -> Result<f32, Error>;
    /// Returns the nearest depth and depth gradient for the given x, y coordinates
    fn depth_and_gradient(&self, x: &f32, y: &f32) -> Result<(f32, (f32, f32)), Error>;
}
