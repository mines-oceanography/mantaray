//! Struct for creating and accessing bathymetry data from an array.
//!
//! This is only used in testing purposes when we purposely want to access out
//! of bounds.

use super::BathymetryData;
use crate::{datatype::{Gradient, Point}, error::Result};

pub(crate) struct ArrayDepth {
    array: Vec<Vec<f32>>,
}

// TODO: to make this `ArrayDepth` useful for use outside generating out of
// bounds values in tests, we need to define grid spacing in both x and y
// directions and map those to cell indexes in the array. Then implement an
// interpolation, and return a valid gradient.
impl BathymetryData for ArrayDepth {
    fn depth(&self, point: &Point<f32>) -> Result<f32> {
        let x = *point.x() as usize;
        let y = *point.y() as usize;
        if x >= self.array.len() || y >= self.array.len() {
            return Ok(f32::NAN);
        }
        Ok(self.array[x][y])
    }

    fn depth_and_gradient(&self, point: &Point<f32>) -> Result<(f32, Gradient<f32>)> {
        let x = *point.x() as usize;
        let y = *point.y() as usize;
        if x >= self.array.len() || y >= self.array.len() {
            return Ok((f32::NAN, Gradient::new(f32::NAN, f32::NAN)));
        }
        Ok((self.array[x][y], Gradient::new(0.0, 0.0)))
    }
}

#[allow(dead_code)]
impl ArrayDepth {
    pub(crate) fn new(array: Vec<Vec<f32>>) -> Self {
        ArrayDepth { array }
    }
}
