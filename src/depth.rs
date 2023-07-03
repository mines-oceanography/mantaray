//! Traits and structs for bathymetry/depth
//! 
//! Contains a trait Depth and two structs ConstantDepth and ArrayDepth which
//! implement the Depth. In the future, ArrayDepth will be changed from a vector
//! of vectors to a more specific type related to the bathymetry data.

use crate::WaveRayPath;

/// A trait with one method to calculate and return the depth for different types.
pub(crate) trait Depth {
    /// Returns the calculated depth at the given x and y.
    fn calculate(&self, x: f64, y: f64) -> f64;
}

impl<'a> WaveRayPath<'a> {
   /// Returns the depth at given x and y.
   pub fn depth(&self, x: f64, y: f64) -> f64 {
      self.data.calculate(x, y)
   }
}

/// A struct that stores the depth as a float and implements the Depth trait.
pub(crate) struct ConstantDepth {
   pub(crate) d: f64,
}

impl Depth for ConstantDepth {
   fn calculate(&self, _x: f64, _y: f64) -> f64 {
       self.d
   }
}

/// A struct that stores the depth as a 2D vector of float and implements the
/// Depth trait.
pub(crate) struct ArrayDepth {
   pub(crate) array: Vec<Vec<f64>>,
}

// FIXME: the program is not crashing, but NAN isn't that useful. maybe give a
// warning for index out of bounds and then set it to NAN
impl Depth for ArrayDepth {
   fn calculate(&self, x: f64, y: f64) -> f64 {
      if x as usize >= self.array.len() || y as usize >= self.array.len() {
         return f64::NAN;
      }
       self.array[x as usize][y as usize] // FIXME: since x and y are floats, they are truncated or rounded to a usize. I probably want a better interpolation estimate
   }
}
