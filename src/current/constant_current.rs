use crate::error::Error;

use super::CurrentData;

pub(crate) struct ConstantCurrent {
    current: f32,
}

impl ConstantCurrent {
    /// Constructor
    pub(crate) fn new(current: f32) -> Self {
        ConstantCurrent { current }
    }
}

impl CurrentData for ConstantCurrent {
    fn get_current(&self, _x: &f32, _y: &f32) -> Result<f32, Error> {
        Ok(self.current)
    }
}
