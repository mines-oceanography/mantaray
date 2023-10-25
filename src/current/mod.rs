use crate::error::Error;

mod constant_current;

#[allow(unused_imports)]
pub(super) use constant_current::ConstantCurrent;

pub(crate) trait CurrentData {
    fn get_current(&self, x: &f32, y: &f32) -> Result<f32, Error>;
}
