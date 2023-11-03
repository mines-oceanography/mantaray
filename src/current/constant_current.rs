use crate::error::Error;

use super::CurrentData;

pub(crate) struct ConstantCurrent {
    u: f32,
    v: f32,
}

#[allow(dead_code)]
impl ConstantCurrent {
    /// Constructor
    pub(crate) fn new(u: f32, v: f32) -> Self {
        ConstantCurrent { u, v }
    }
}

impl CurrentData for ConstantCurrent {
    fn current(&self, _x: &f32, _y: &f32) -> Result<(f32, f32), Error> {
        Ok((self.u, self.v))
    }
}
