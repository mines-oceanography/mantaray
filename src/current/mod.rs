//! CurrentData
//!
//! This module contains the following structs that implement the `CurrentData`
//! trait:
//! - `ConstantCurrent`

use crate::datatype::{Current, Gradient, Point};
use crate::error::Result;

mod cartesian_current;
mod constant_current;

#[allow(unused_imports)]
pub(super) use cartesian_current::CartesianCurrent;
#[allow(unused_imports)]
pub(super) use constant_current::ConstantCurrent;
#[allow(unused_imports)]
pub(super) use constant_current::DEFAULT_CURRENT;

/// A trait implementing methods to get current and gradient
pub(crate) trait CurrentData: Sync {
    #[allow(dead_code)]
    /// Current (u, v) at the given (x, y)
    fn current(&self, point: &Point<f64>) -> Result<Current<f64>>;

    /// Current (u, v) and the gradient (du/dx, du/dy, dv/dx, dv/dy)
    fn current_and_gradient(
        &self,
        point: &Point<f64>,
    ) -> Result<(Current<f64>, (Gradient<f64>, Gradient<f64>))>;
}
