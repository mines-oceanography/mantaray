use crate::error::Error;

use super::CurrentData;

pub(crate) struct FunctionCurrent {
    /// `fn_current` is a function pointer to a function with the signature
    /// `fn(&f64, &f64) -> Result<(f64, f64), Error>`. This function represents
    /// taking an input of (x, y) and returning either a (u, v) or an error.
    fn_current: fn(&f64, &f64) -> Result<(f64, f64), Error>,
}

#[allow(dead_code)]
impl FunctionCurrent {
    /// create a new instance of `FunctionCurrent` with a given function
    /// pointer.
    ///
    /// # Arguments:
    ///
    /// - `fn_current`: `fn(&f64, &f64) -> Result<(f64, f64), Error>`
    ///
    /// The function that calculates current given an (x, y) input. A function
    /// that takes two references to f64 values and returns a Result containing
    /// a tuple of two f64 values or an Error.
    ///
    /// # Returns:
    ///
    /// `FunctionCurrent` : returns the constructed struct
    pub(crate) fn new(fn_current: fn(&f64, &f64) -> Result<(f64, f64), Error>) -> Self {
        FunctionCurrent { fn_current }
    }
}

impl CurrentData for FunctionCurrent {
    /// Return (u, v) at given (x, y)
    ///
    /// # Arguments
    /// - `x` : `&f64` the x coordinate
    ///
    /// - `y` : `&f64` the y coordinate
    ///
    /// # Return
    /// `Result<(f64, f64), Error>` : return either (u, v) or an error.
    ///
    /// # Error
    /// `current` can return an error that appears in the function `fn_current`
    /// references.
    ///
    /// # Panics
    /// `current` can panic if the function in `fn_current` panics.
    fn current(&self, x: &f64, y: &f64) -> Result<(f64, f64), Error> {
        (self.fn_current)(x, y)
    }

    /// Return (u, v), (du/dx, du/dy, dv/dx, dv/dy) at given (x, y)
    ///
    /// # Arguments
    /// - `x` : `&f64` the x coordinate
    ///
    /// - `y` : `&f64` the y coordinate
    ///
    /// # Return
    /// `Result<((f64, f64), (f64, f64, f64, f64)), Error>` : return either ((u,
    /// v), (du/dx, du/dy, dv/dx, dv/dy)) or an error.
    ///
    /// # Error
    /// `current` can return an error that appears in the function `fn_current`
    /// references.
    ///
    /// # Panics
    /// `current` can panic if the function in `fn_current` panics.
    fn current_and_gradient(
        &self,
        x: &f64,
        y: &f64,
    ) -> Result<((f64, f64), (f64, f64, f64, f64)), Error> {
        let (u, v) = (self.fn_current)(x, y)?;

        let d = 1.0; // FIXME: this needs to be a better estimate?
        let (u_dx, v_dx) = self.current(&(x + d), y)?;
        let (u_dy, v_dy) = self.current(x, &(y + d))?;
        let (u_dx_i, v_dx_i) = self.current(&(x - d), y)?;
        let (u_dy_i, v_dy_i) = self.current(x, &(y - d))?;

        let du_dx = (u_dx - u_dx_i) / (2.0 * d);
        let du_dy = (u_dy - u_dy_i) / (2.0 * d);
        let dv_dx = (v_dx - v_dx_i) / (2.0 * d);
        let dv_dy = (v_dy - v_dy_i) / (2.0 * d);

        Ok(((u, v), (du_dx, du_dy, dv_dx, dv_dy)))
    }
}

#[cfg(test)]
mod test_constant_change_current {
    use crate::{
        current::{CurrentData, FunctionCurrent},
        error::Error,
    };

    #[test]
    fn test_current_x() {
        fn fn_current(x: &f64, _y: &f64) -> Result<(f64, f64), Error> {
            let u = 1.0 + x;
            let v = 1.0 + x;
            Ok((u, v))
        }

        let current_data = FunctionCurrent::new(fn_current);

        // simple test case x
        // check at at (0,0) returns (u0, v0)
        assert!((current_data.current(&0.0, &0.0).unwrap().0 - 1.0).abs() < f64::EPSILON);
        assert!((current_data.current(&0.0, &0.0).unwrap().1 - 1.0).abs() < f64::EPSILON);
        // check at (1, 1) returns (u0 + 1, v0 + 1)
        assert!((current_data.current(&1.0, &1.0).unwrap().0 - 2.0).abs() < f64::EPSILON);
        assert!((current_data.current(&1.0, &1.0).unwrap().1 - 2.0).abs() < f64::EPSILON);
        // check at (3, 5) returns (u0 + 3, v0 + 3)
        assert!((current_data.current(&3.0, &5.0).unwrap().0 - 4.0).abs() < f64::EPSILON);
        assert!((current_data.current(&3.0, &5.0).unwrap().1 - 4.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_current_y() {
        fn fn_current(_x: &f64, y: &f64) -> Result<(f64, f64), Error> {
            let u = 1.0 + y;
            let v = 1.0 + y;
            Ok((u, v))
        }

        let current_data = FunctionCurrent::new(fn_current);

        // simple test case x
        // check at at (0,0) returns (u0, v0)
        assert!((current_data.current(&0.0, &0.0).unwrap().0 - 1.0).abs() < f64::EPSILON);
        assert!((current_data.current(&0.0, &0.0).unwrap().1 - 1.0).abs() < f64::EPSILON);
        // check at (1, 1) returns (u0 + 1, v0 + 1)
        assert!((current_data.current(&1.0, &1.0).unwrap().0 - 2.0).abs() < f64::EPSILON);
        assert!((current_data.current(&1.0, &1.0).unwrap().1 - 2.0).abs() < f64::EPSILON);
        // check at (3, 5) returns (u0 + 3, v0 + 3)
        assert!((current_data.current(&3.0, &5.0).unwrap().0 - 6.0).abs() < f64::EPSILON);
        assert!((current_data.current(&3.0, &5.0).unwrap().1 - 6.0).abs() < f64::EPSILON);
    }

    #[test]
    /// FIXME: this is going to be an estimate, so we need to determine the best
    /// way to estimate the gradient of u and v.
    fn test_gradient() {
        fn fn_current(x: &f64, _y: &f64) -> Result<(f64, f64), Error> {
            let u = 1.0 + x;
            let v = 1.0 + x;
            Ok((u, v))
        }

        let current_data = FunctionCurrent::new(fn_current);

        // simple test case x
        // check at at (0,0) returns (u0, v0)
        let ((u, v), (du_dx, du_dy, dv_dx, dv_dy)) =
            current_data.current_and_gradient(&0.0, &0.0).unwrap();
        assert!((u - 1.0).abs() < f64::EPSILON);
        assert!((v - 1.0).abs() < f64::EPSILON);
        assert!((du_dx - 1.0).abs() < f64::EPSILON);
        assert!((du_dy - 0.0).abs() < f64::EPSILON);
        assert!((dv_dx - 1.0).abs() < f64::EPSILON);
        assert!((dv_dy - 0.0).abs() < f64::EPSILON);

        // check at (1, 1) returns (u0 + 1, v0 + 1)
        let ((u, v), (du_dx, du_dy, dv_dx, dv_dy)) =
            current_data.current_and_gradient(&1.0, &1.0).unwrap();
        assert!((u - 2.0).abs() < f64::EPSILON);
        assert!((v - 2.0).abs() < f64::EPSILON);
        assert!((du_dx - 1.0).abs() < f64::EPSILON);
        assert!((du_dy - 0.0).abs() < f64::EPSILON);
        assert!((dv_dx - 1.0).abs() < f64::EPSILON);
        assert!((dv_dy - 0.0).abs() < f64::EPSILON);
    }
}
