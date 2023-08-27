use super::BathymetryData;
use crate::error::Error;
use derive_builder::Builder;

#[derive(Builder, Debug, PartialEq)]
pub(crate) struct ConstantDepth {
    #[builder(default = "1000.0")]
    h: f32,
}

impl BathymetryData for ConstantDepth {
    fn get_depth(&self, _x: &f32, _y: &f32) -> Result<f32, Error> {
        Ok(self.h)
    }

    fn get_depth_and_gradient(&self, _x: &f32, _y: &f32) -> Result<(f32, (f32, f32)), Error> {
        Ok((self.h, (0.0, 0.0)))
    }
}

impl ConstantDepth {
    fn builder() -> ConstantDepthBuilder {
        ConstantDepthBuilder::default()
    }

    pub(crate) fn new(h: f32) -> ConstantDepth {
        ConstantDepth { h }
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
