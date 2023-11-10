//! CurrentData
//! 
//! This module contains the following structs that implement the `CurrentData`
//! trait:
//! - `ConstantCurrent`

use crate::error::Error;

mod cartesian_file_current;
mod constant_current;

#[allow(unused_imports)]
pub(super) use constant_current::ConstantCurrent;

pub(crate) trait CurrentData {
    /// Return the current (u, v) at the given (x, y)
    fn current(&self, x: &f64, y: &f64) -> Result<(f64, f64), Error>;
}
