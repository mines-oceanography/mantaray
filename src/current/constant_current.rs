use crate::error::Error;

use super::CurrentData;

pub(crate) struct ConstantCurrent {
    u: f64,
    v: f64,
}

#[allow(dead_code)]
impl ConstantCurrent {
    /// Constructor
    pub(crate) fn new(u: f64, v: f64) -> Self {
        ConstantCurrent { u, v }
    }
}

impl CurrentData for ConstantCurrent {
    fn current(&self, _x: &f64, _y: &f64) -> Result<(f64, f64), Error> {
        Ok((self.u, self.v))
    }
}
