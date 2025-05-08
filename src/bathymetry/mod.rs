//! Trait and structs for accessing depth and gradient from various bathymetry
//! data.
//!
//! The implementors of the `BathymetryData` trait are different types of
//! bathymetry:
//! - `CartesianNetcdf3` - read and access the data stored in a NetCDF3 file.
//! - `ConstantDepth` - constant depth bathymetry. There are no domain
//!   constraints on the input since the depth is defined by a constant value.
//! - `ConstantSlope` - constant slope bathymetry. There are no domain
//!   constraints on the input since the depth is defined by a function.
//!
//! The following are used primarily for testing purposes:
//! - `ArrayDepth` - used to create bathymetry data from an array. Useful for
//!   creating purposefully out of bounds points.

mod array_depth;
mod cartesian_netcdf3;
mod constant_depth;
mod constant_slope;

use crate::datatype::{Gradient, Point};
use crate::error::Result;
#[allow(unused_imports)]
pub(super) use array_depth::ArrayDepth;
#[allow(unused_imports)]
pub(super) use cartesian_netcdf3::CartesianNetcdf3;
#[allow(unused_imports)]
pub(super) use constant_depth::ConstantDepth;
#[allow(unused_imports)]
pub(super) use constant_depth::DEFAULT_BATHYMETRY;
#[allow(unused_imports)]
pub(super) use constant_slope::ConstantSlope;

/// A trait defining ability to return depth and gradient
pub(crate) trait BathymetryData: Sync {
    #[allow(dead_code)]
    /// Returns the nearest depth for the given (x, y) point.
    fn depth(&self, point: &Point<f32>) -> Result<f32>;
    /// Returns the nearest depth and depth gradient for the given (x, y) coordinates
    fn depth_and_gradient(&self, point: &Point<f32>) -> Result<(f32, Gradient<f32>)>;
}
