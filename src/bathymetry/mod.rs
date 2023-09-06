//! Bathymetry

pub(crate) mod cartesian;
mod constant_depth;
mod constant_slope;

use crate::error::Error;
#[allow(unused_imports)]
pub(super) use cartesian::CartesianFile;
#[allow(unused_imports)]
pub(super) use constant_depth::ConstantDepth;
#[allow(unused_imports)]
pub(super) use constant_slope::ConstantSlope;

/// A trait used to give the function get_depth
pub(crate) trait BathymetryData: Sync {
    /// Returns the nearest depth for the given x, y coordinate.
    fn get_depth(&self, x: &f32, y: &f32) -> Result<f32, Error>;
    /// Returns the nearest depth and depth gradient for the given x, y coordinates
    fn get_depth_and_gradient(&self, x: &f32, y: &f32) -> Result<(f32, (f32, f32)), Error>;
}
