use crate::error::Error;

use super::CurrentData;
use derive_builder::Builder;

#[derive(Builder, Debug, PartialEq)]
/// The `ConstantChange` struct represents current field with a constant rate of
/// change. The user sets initial values of (x, y), (u, v), as well as a
/// constant value for the gradients of u and v.
///
/// # Properties:
///
///`x0`: `f64`
/// - Initial value of x
///
///`u0`: `f64`
/// - The `u0` property represents the initial value of u.
///
///`dudx`: `f64`
/// - The property `dudx` represents the rate of change of `u` with respect to
///   `x`.
///
///`dudy`: `f64`
/// - The property `dudy` represents the rate of change of `u` with respect to
///   `y`.
///
///`y0`: `f64`
/// - The property `y0` represents the initial value of the variable `y`.
///
///`v0`: `f64`
/// - The property `v0` represents the initial value of the v.
///
///`dvdx`: `f64`
/// - The property `dvdx` represents the rate of change of `v0` with respect to
///   `x`.
///
///`dvdy`: `f64`
/// - The `dvdy` property represents the rate of change of `v0` with respect to
///   `y`.
pub(crate) struct ConstantChange {
    #[builder(default = "0.0")]
    /// The initial x value
    x0: f64,
    #[builder(default = "0.0")]
    /// the initial u value
    u0: f64,
    #[builder(default = "0.0")]
    /// the partial of u with respect to x
    dudx: f64,
    #[builder(default = "0.0")]
    /// the partial of v with respect to y
    dudy: f64,
    #[builder(default = "0.0")]
    /// the initial y value
    y0: f64,
    #[builder(default = "0.0")]
    /// the initial v value
    v0: f64,
    #[builder(default = "0.0")]
    /// the partial of v with respect to x
    dvdx: f64,
    #[builder(default = "0.0")]
    /// the partial of v with respect to64
    dvdy: f64,
}

#[allow(dead_code)]
impl ConstantChange {
    /// Constructor creates a new instance of the `ConstantSlope` struct with
    /// the given initial values.
    ///
    /// Arguments:
    ///
    /// * `x0`: `f64`
    /// The initial x-coordinate.
    ///
    /// * `y0`: `f64`
    /// The initial y-coordinate.
    ///
    /// * `u0`: `f64`
    /// Initial value of the x-component of the current vector.
    ///
    /// * `dudx`: `f64`
    /// The parameter `dudx` represents the rate of change of the
    /// x-component of the current with respect to x.
    ///
    /// * `dudy`: `f64`
    /// The parameter `dudy` represents the rate of change of the
    /// x-component of the current with respect to y.
    ///
    ///  * `v0`: `f64`
    /// The parameter `v0` represents the initial value of the
    /// y-component of the current.
    ///
    /// * `dvdx`: `f64`
    /// The parameter `dvdx` represents the rate of change of the
    ///  y-component of the current with respect to x.
    ///
    /// * `dvdy`: `f64`
    /// The parameter `dvdy` represents the rate of change of the
    /// y-component of the current with respect to changes in y.
    ///
    /// Returns:
    ///
    /// The constructed `ConstantSlope` struct is returned.
    pub(crate) fn new(
        x0: f64,
        y0: f64,
        u0: f64,
        dudx: f64,
        dudy: f64,
        v0: f64,
        dvdx: f64,
        dvdy: f64,
    ) -> Self {
        ConstantChange {
            u0,
            dudx,
            dudy,
            v0,
            dvdx,
            dvdy,
            x0,
            y0,
        }
    }

    /// `build` returns a `ConstantChangeBuilder` object.
    ///
    /// # Returns:
    /// - `ConstantChangeBuilder::default()`
    pub(crate) fn build() -> ConstantChangeBuilder {
        ConstantChangeBuilder::default()
    }
}

impl CurrentData for ConstantChange {
    // Return (u, v)
    fn current(&self, x: &f64, y: &f64) -> Result<(f64, f64), Error> {
        let u = self.u0 + (x - self.x0) * self.dudx + (y - self.y0) * self.dudy;
        let v = self.v0 + (x - self.x0) * self.dvdx + (y - self.y0) * self.dvdy;
        Ok((u, v))
    }

    // Return (u, v), (du/dx, du/dy, dv/dx, dv/dy)
    fn current_and_gradient(
        &self,
        x: &f64,
        y: &f64,
    ) -> Result<((f64, f64), (f64, f64, f64, f64)), Error> {
        let u = self.u0 + (x - self.x0) * self.dudx + (y - self.y0) * self.dudy;
        let v = self.v0 + (x - self.x0) * self.dvdx + (y - self.y0) * self.dvdy;
        Ok(((u, v), (self.dudx, self.dudy, self.dvdx, self.dvdy)))
    }
}

mod test_constant_change_current {
    use crate::current::{ConstantChange, CurrentData};

    #[test]
    fn test_current() {
        // simple test case x
        // x0, y0, u0, dudx, dudy, v0, dvdx, dvdy
        let (x0, y0, u0, dudx, dudy, v0, dvdx, dvdy) = (0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0);
        let current_data = ConstantChange::new(x0, y0, u0, dudx, dudy, v0, dvdx, dvdy);
        // check at at (0,0) returns (u0, v0)
        assert!((current_data.current(&0.0, &0.0).unwrap().0 - u0).abs() < f64::EPSILON);
        assert!((current_data.current(&0.0, &0.0).unwrap().1 - v0).abs() < f64::EPSILON);
        // check at (1, 1) returns (u0 + 1, v0 + 1)
        assert!((current_data.current(&1.0, &1.0).unwrap().0 - 2.0).abs() < f64::EPSILON);
        assert!((current_data.current(&1.0, &1.0).unwrap().1 - 2.0).abs() < f64::EPSILON);
        // check at (3, 5) returns (u0 + 3, v0 + 3)
        assert!((current_data.current(&3.0, &5.0).unwrap().0 - 4.0).abs() < f64::EPSILON);
        assert!((current_data.current(&3.0, &5.0).unwrap().1 - 4.0).abs() < f64::EPSILON);

        // simple test case y
        let (x0, y0, u0, dudx, dudy, v0, dvdx, dvdy) = (0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0);
        let current_data = ConstantChange::new(x0, y0, u0, dudx, dudy, v0, dvdx, dvdy);
        // check at at (0,0) returns (u0, v0)
        assert!((current_data.current(&0.0, &0.0).unwrap().0 - 1.0).abs() < f64::EPSILON);
        assert!((current_data.current(&0.0, &0.0).unwrap().1 - 1.0).abs() < f64::EPSILON);
        // check at (1, 1) returns (u0 + 1, v0 + 1)
        assert!((current_data.current(&1.0, &1.0).unwrap().0 - 2.0).abs() < f64::EPSILON);
        assert!((current_data.current(&1.0, &1.0).unwrap().1 - 2.0).abs() < f64::EPSILON);
        // check at (3, 5) returns (u0 + 5, v0 + 5)
        assert!((current_data.current(&3.0, &5.0).unwrap().0 - 6.0).abs() < f64::EPSILON);
        assert!((current_data.current(&3.0, &5.0).unwrap().1 - 6.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_builder() {
        let current_data = ConstantChange::build().dudx(0.0).dvdx(0.0).build().unwrap();
        assert_eq!(
            current_data,
            ConstantChange::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        )
    }
}
