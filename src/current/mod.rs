use crate::error::Error;

mod constant_current;

#[allow(unused_imports)]
pub(super) use constant_current::ConstantCurrent;

pub(crate) trait CurrentData {
    /// Return the current (u, v) at the given (x, y)
    fn current(&self, x: &f32, y: &f32) -> Result<(f32, f32), Error>;
}
