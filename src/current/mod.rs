//! CurrentData
//!
//! This module contains the following structs that implement the `CurrentData`
//! trait:
//! - `ConstantCurrent`

use crate::error::Result;
use crate::{Current, Point};

mod cartesian_current;
mod constant_current;

#[allow(unused_imports)]
pub(super) use cartesian_current::CartesianCurrent;
#[allow(unused_imports)]
pub(super) use constant_current::ConstantCurrent;

pub trait CurrentData: Sync {
    /// Current (u, v) at the given (x, y)
    fn current(&self, point: &Point<f64>) -> Result<Current<f64>>;

    /// Current (u, v) and the gradient (du/dx, du/dy, dv/dx, dv/dy)
    fn current_and_gradient(
        &self,
        point: &Point<f64>,
    ) -> Result<(Current<f64>, (f64, f64, f64, f64))>;
}
