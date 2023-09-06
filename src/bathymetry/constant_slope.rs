use super::BathymetryData;
use crate::error::Error;
use derive_builder::Builder;

#[derive(Builder, Debug, PartialEq)]
/// A bathymetry database with constant slope
///
/// This might be only useful for development and tests.
pub(crate) struct ConstantSlope {
    #[builder(default = "50.0")]
    h0: f32,
    #[builder(default = "0.0")]
    x0: f32,
    #[builder(default = "0.0")]
    y0: f32,
    #[builder(default = "5e-2")]
    dx: f32,
    #[builder(default = "0.0")]
    dy: f32,
}

impl BathymetryData for ConstantSlope {
    /// Depth for a given position (x, y)
    ///
    /// Returns NaN when any input is NaN. Since it is a constant slope,
    /// there is no concept of boundaries, thus it can't fail as out of
    /// bounds.
    fn get_depth(&self, x: &f32, y: &f32) -> Result<f32, Error> {
        if x.is_nan() || y.is_nan() {
            Ok(f32::NAN)
        } else {
            Ok(self.h0 + self.dx * (x - self.x0) + self.dy * (y - self.y0))
        }
    }

    /// Depth and gradient for a given position (x, y)
    ///
    /// Returns NaN when any input is NaN. Since it is a constant slope,
    /// there is no concept of boundaries, thus it can't fail as out of
    /// bounds.
    fn get_depth_and_gradient(&self, x: &f32, y: &f32) -> Result<(f32, (f32, f32)), Error> {
        if x.is_nan() || y.is_nan() {
            Ok((f32::NAN, (f32::NAN, f32::NAN)))
        } else {
            let h = self.h0 + self.dx * (x - self.x0) + self.dy * (y - self.y0);
            Ok((h, (self.dx, self.dy)))
        }
    }
}

impl ConstantSlope {
    #[allow(dead_code)]
    fn builder() -> ConstantSlopeBuilder {
        ConstantSlopeBuilder::default()
    }
}

#[cfg(test)]
mod test_constant_slope {
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
            dx: 1e-2,
            dy: 0.0,
        };

        assert!(c.get_depth(&f32::NAN, &0.0).unwrap().is_nan());
        assert!(c.get_depth(&0.0, &f32::NAN).unwrap().is_nan());
        assert!(c.get_depth(&f32::NAN, &f32::NAN).unwrap().is_nan());
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
                h0: 1000.0,
                x0: 0.0,
                y0: 0.0,
                dx: 1e-2,
                dy: 0.0
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
                dx: 1e-2,
                dy: 0.0
            }
        );
    }

    #[test]
    fn builder() {
        let c = ConstantSlope::builder().build().unwrap();
        assert_eq!(
            c,
            ConstantSlope {
                h0: 1000.0,
                x0: 0.0,
                y0: 0.0,
                dx: 1e-2,
                dy: 0.0
            }
        );
    }
}
