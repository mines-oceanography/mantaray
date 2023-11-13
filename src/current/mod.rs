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

pub(crate) trait CurrentData {
    /// Return the current (u, v) at the given (x, y)
    fn current(&self, x: &f64, y: &f64) -> Result<(f64, f64), Error>;
    fn current_and_gradient(
        &self,
        x: &f64,
        y: &f64,
    ) -> Result<((f64, f64), (f64, f64, f64, f64)), Error>;
}
