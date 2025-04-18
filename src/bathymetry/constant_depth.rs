//! Struct used to create and access bathymetry data with a constant depth.

use super::BathymetryData;
use crate::{datatype::{Gradient, Point}, error::Result};
use derive_builder::Builder;

#[allow(dead_code)]
/// The default depth is 2000 m
pub(crate) const DEFAULT_BATHYMETRY: ConstantDepth = ConstantDepth { h: 2000.0 };

#[derive(Builder, Debug, PartialEq)]
/// A bathymetry database with constant depth
///
/// This might be only useful for development and tests.
pub(crate) struct ConstantDepth {
    #[builder(default = "1000.0")]
    h: f32,
}

impl BathymetryData for ConstantDepth {
    /// Depth for a given position (x, y)
    ///
    /// Returns NaN when any input is NaN. Since it is a constant depth,
    /// there is no concept of boundaries, thus it can't fail as out of
    /// bounds.
    fn depth(&self, point: &Point<f32>) -> Result<f32> {
        if point.x().is_nan() || point.y().is_nan() {
            Ok(f32::NAN)
        } else {
            Ok(self.h)
        }
    }

    /// Depth and gradient for a given position (x, y)
    ///
    /// Returns NaN when any input is NaN. Since it is a constant depth,
    /// there is no concept of boundaries, thus it can't fail as out of
    /// bounds.
    fn depth_and_gradient(&self, point: &Point<f32>) -> Result<(f32, Gradient<f32>)> {
        if point.x().is_nan() || point.y().is_nan() {
            Ok((f32::NAN, Gradient::new(f32::NAN, f32::NAN)))
        } else {
            Ok((self.h, Gradient::new(0.0, 0.0)))
        }
    }
}

impl ConstantDepth {
    #[allow(dead_code)]
    fn builder() -> ConstantDepthBuilder {
        ConstantDepthBuilder::default()
    }

    #[allow(dead_code)]
    /// construct a bathymetry dataset with constant depth
    ///
    /// # Arguments
    /// `h`: `f32` the depth [m]
    pub(crate) fn new(h: f32) -> ConstantDepth {
        ConstantDepth { h }
    }
}

#[cfg(test)]
mod test_constant_depth {
    use crate::datatype::Point;

    use super::{BathymetryData, ConstantDepth};

    #[test]
    // Maybe this is not clear for constant depth, but it should respect the
    // general behavior for other types of bathymetries. A NaN input is not
    // an error per se, but should result in a NaN result. Different than
    // an out of bounds request or unfeasible lat/lon coordinate.
    fn nan_input() {
        let c = ConstantDepth { h: 100.0 };

        assert!(c.depth(&Point::new(f32::NAN, 0.0)).unwrap().is_nan());
        assert!(c.depth(&Point::new(0.0, f32::NAN)).unwrap().is_nan());
        assert!(c.depth(&Point::new(f32::NAN, f32::NAN)).unwrap().is_nan());
    }
}

#[cfg(test)]
mod test_builder {
    use super::{ConstantDepth, ConstantDepthBuilder};

    #[test]
    fn build_default() {
        let c = ConstantDepthBuilder::default().build().unwrap();
        assert_eq!(c, ConstantDepth { h: 1000.0 });
    }

    #[test]
    fn build() {
        let c = ConstantDepthBuilder::default().h(42.0).build().unwrap();
        assert_eq!(c, ConstantDepth { h: 42.0 });
    }

    #[test]
    fn builder() {
        let c = ConstantDepth::builder().build().unwrap();
        assert_eq!(c, ConstantDepth { h: 1000.0 });
    }
}
