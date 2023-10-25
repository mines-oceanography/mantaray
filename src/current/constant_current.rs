use crate::error::Error;

use super::CurrentData;

pub(crate) struct ConstantCurrent {
    ux: f32,
    uy: f32,
}

impl ConstantCurrent {
    /// Constructor
    pub(crate) fn new(ux: f32, uy: f32) -> Self {
        ConstantCurrent { ux, uy }
    }
}

impl CurrentData for ConstantCurrent {
    fn current(&self, _x: &f32, _y: &f32) -> Result<(f32, f32), Error> {
        Ok((self.ux, self.uy))
    }
}
