use crate::error::Result;

use super::CurrentData;

/// struct representing a constant current field.
///
/// # Properties:
///
/// * `u`: `f64` value representing x component of the current.
///
/// * `v`: `f64` value representing y component of the current.
pub(crate) struct ConstantCurrent {
    /// x component of the current
    u: f64,
    /// y component of the current
    v: f64,
}

#[allow(dead_code)]
impl ConstantCurrent {
    /// Constructor
    ///
    /// # Returns
    /// returns the constructed ConstantCurrent
    pub(crate) fn new(u: f64, v: f64) -> Self {
        ConstantCurrent { u, v }
    }
}

impl CurrentData for ConstantCurrent {
    /// get the current at point (x, y)
    ///
    /// # Arguments
    /// - `x` : `f64` the x location
    ///
    /// - `y` : `f64` the y location
    ///
    /// # Returns
    /// `Result<(f64, f64), Error>` : returns the values (u, v) or an Error.
    ///
    /// # Error
    /// The trait definition includes the chance for error. However, the
    /// `ConstantCurrent::current` should never return an error.
    fn current(&self, _x: &f64, _y: &f64) -> Result<(f64, f64)> {
        Ok((self.u, self.v))
    }

    /// get the current and gradient at point (x, y)
    ///
    /// # Arguments
    /// - `x` : `f64` the x location
    ///
    /// - `y` : `f64` the y location
    ///
    /// # Returns
    /// `Result<((f64, f64), (f64, f64, f64, f64)), Error>` : returns the values
    /// (u, v) and (du/dx, du/dy, dv/dx, dv/dy) or an Error.
    ///
    /// # Error
    /// The trait definition includes the chance for error. However, the
    /// `ConstantCurrent::current_and_gradient` should never return an error.
    fn current_and_gradient(
        &self,
        _x: &f64,
        _y: &f64,
    ) -> Result<((f64, f64), (f64, f64, f64, f64))> {
        Ok(((self.u, self.v), (0.0, 0.0, 0.0, 0.0)))
    }
}
