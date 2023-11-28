//! CurrentData
//!
//! This module contains the following structs that implement the `CurrentData`
//! trait:
//! - `ConstantCurrent`

use crate::error::Error;

mod constant_current;
mod constant_rate_change_current;
mod function_current;

#[allow(unused_imports)]
pub(super) use constant_current::ConstantCurrent;
#[allow(unused_imports)]
pub(super) use constant_rate_change_current::ConstantChange;
#[allow(unused_imports)]
pub(super) use function_current::FunctionCurrent;

#[allow(dead_code)]
/// The CurrentType is a type alias for a tuple of two f64 values. This tuple
/// represents (u, v) current.
pub(crate) type CurrentType = (f64, f64);
#[allow(dead_code)]
/// The GradientType is a type alias for a tuple of four f64 values. This tuple
/// represents (du/dx, du/dy, dv/dx, dv/dy), or the gradients of u and v.
pub(crate) type GradientType = (f64, f64, f64, f64);

pub(crate) trait CurrentData {
    /// Return the current (u, v) at the given (x, y)
    fn current(&self, x: &f64, y: &f64) -> Result<(f64, f64), Error>;
    fn current_and_gradient(&self, x: &f64, y: &f64) -> Result<(CurrentType, GradientType), Error>;
}
