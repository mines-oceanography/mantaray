//! CurrentData
//!
//! This module contains the following structs that implement the `CurrentData`
//! trait:
//! - `ConstantCurrent`

use crate::error::Result;

mod cartesian_current;
mod constant_current;

#[allow(unused_imports)]
pub(super) use cartesian_current::CartesianCurrent;
#[allow(unused_imports)]
pub(super) use constant_current::ConstantCurrent;

pub trait CurrentData: Sync {
    /// Return the current (u, v) at the given (x, y)
    fn current(&self, x: &f64, y: &f64) -> Result<(f64, f64)>;
    /// Return the current (u, v) and the gradient (du/dx, du/dy, dv/dx, dv/dy)
    fn current_and_gradient(&self, x: &f64, y: &f64) -> Result<((f64, f64), (f64, f64, f64, f64))>;
}
