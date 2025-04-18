//! Struct used to create and access bathymetry data with a constant slope.

use super::BathymetryData;
use crate::{datatype::{Gradient, Point}, error::Result};
use derive_builder::Builder;

#[derive(Builder, Debug, PartialEq)]
/// A bathymetry pseudo-database with constant slope
///
/// # Arguments
///
/// * `h0`: Depth [m] at ($x_0$, $y_0$).
/// * `x0`: `x` coordinate [m] where depth is `h0` at `y0`.
/// * `y0`: `y` coordinate [m] where depth is `h0` at `x0`.
/// * `dhdx`: Slope on x direction.
/// * `dhdy`: Slope on y direction.
///
/// ConstantSlope can be build with default values, which is a convenient
/// approach to build tests.
///
/// * Manual definition
/// let bathymetry = ConstantSlope(100, 0, 0, -1e-2, 0);
///
/// * Using default values
/// let h100 = ConstantSlope::builder().h0(100).build().unwrap();
///
/// This might be only useful for development and tests.
pub(crate) struct ConstantSlope {
    /// initial depth
    #[builder(default = "50.0")]
    h0: f32,
    ///initial x
    #[builder(default = "0.0")]
    x0: f32,
    ///initial y
    #[builder(default = "0.0")]
    y0: f32,
    ///rate of change in depth with respect to x
    #[builder(default = "-5e-2")]
    dhdx: f32,
    ///rate of change in depth with respect to y
    #[builder(default = "0.0")]
    dhdy: f32,
}

impl BathymetryData for ConstantSlope {
    /// Depth for a given position (x, y)
    ///
    /// Returns NaN when any input is NaN. Since it is a constant slope,
    /// there is no concept of boundaries, thus it can't fail as out of
    /// bounds.
    fn depth(&self, point: &Point<f32>) -> Result<f32> {
        let x = point.x();
        let y = point.y();
        if x.is_nan() || y.is_nan() {
            Ok(f32::NAN)
        } else {
            Ok(self.h0 + self.dhdx * (x - self.x0) + self.dhdy * (y - self.y0))
        }
    }

    /// Depth and gradient for a given position (x, y)
    ///
    /// Returns NaN when any input is NaN. Since it is a constant slope,
    /// there is no concept of boundaries, thus it can't fail as out of
    /// bounds.
    fn depth_and_gradient(&self, point: &Point<f32>) -> Result<(f32, Gradient<f32>)> {
        let x = point.x();
        let y = point.y();
        if x.is_nan() || y.is_nan() {
            Ok((f32::NAN, Gradient::new(f32::NAN, f32::NAN)))
        } else {
            let h = self.h0 + self.dhdx * (x - self.x0) + self.dhdy * (y - self.y0);
            Ok((h, Gradient::new(self.dhdx, self.dhdy)))
        }
    }
}

impl ConstantSlope {
    #[allow(dead_code)]
    /// create the default `ConstantSlopeBuilder` object
    ///
    /// For example,
    /// `ConstantSlope::builder().h0(100.0).dhdx(-0.05).build().unwrap()` builds
    /// the default ConstantSlope, but sets the initial height to 100 m and the
    /// rate of change of depth in the x direction to -0.05.
    pub(crate) fn builder() -> ConstantSlopeBuilder {
        ConstantSlopeBuilder::default()
    }
}

#[cfg(test)]
mod test_constant_slope {
    use crate::datatype::Point;

    use super::{BathymetryData, ConstantSlope};

    #[test]
    // Maybe this is not clear for constant slope, but it should respect the
    // general behavior for other types of bathymetries. A NaN input is not
    // an error per se, but should result in a NaN result. Different than
    // an out of bounds request or unfeasible lat/lon coordinate.
    fn nan_input() {
        let c = ConstantSlope {
            h0: 100.0,
            x0: 0.0,
            y0: 0.0,
            dhdx: -1e-2,
            dhdy: 0.0,
        };

        assert!(c.depth(&Point::new(f32::NAN, 0.0)).unwrap().is_nan());
        assert!(c.depth(&Point::new(0.0, f32::NAN)).unwrap().is_nan());
        assert!(c.depth(&Point::new(f32::NAN, f32::NAN)).unwrap().is_nan());
    }
}

#[cfg(test)]
mod test_builder {
    use super::{ConstantSlope, ConstantSlopeBuilder};

    #[test]
    fn build_default() {
        let c = ConstantSlopeBuilder::default().build().unwrap();
        assert_eq!(
            c,
            ConstantSlope {
                h0: 50.0,
                x0: 0.0,
                y0: 0.0,
                dhdx: -5e-2,
                dhdy: 0.0
            }
        );
    }

    #[test]
    fn build() {
        let c = ConstantSlopeBuilder::default().h0(42.0).build().unwrap();
        assert_eq!(
            c,
            ConstantSlope {
                h0: 42.0,
                x0: 0.0,
                y0: 0.0,
                dhdx: -5e-2,
                dhdy: 0.0
            }
        );
    }

    #[test]
    fn builder() {
        let c = ConstantSlope::builder().build().unwrap();
        assert_eq!(
            c,
            ConstantSlope {
                h0: 50.0,
                x0: 0.0,
                y0: 0.0,
                dhdx: -5e-2,
                dhdy: 0.0
            }
        );
    }
}
