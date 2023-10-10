use super::BathymetryData;
use crate::error::Error;

pub(crate) struct ArrayDepth {
    array: Vec<Vec<f32>>,
}

// FIXME: the program is not crashing, but NAN isn't that useful. maybe use option: some or none
impl BathymetryData for ArrayDepth {
    fn get_depth(&self, x: &f32, y: &f32) -> Result<f32, Error> {
        if *x as usize >= self.array.len() || *y as usize >= self.array.len() {
            return Ok(f32::NAN);
        }
        Ok(self.array[*x as usize][*y as usize]) // FIXME: since x and y are floats, they are truncated or rounded to a usize. I probably want a better interpolation estimate
    }

    fn get_depth_and_gradient(&self, x: &f32, y: &f32) -> Result<(f32, (f32, f32)), Error> {
        if *x as usize >= self.array.len() || *y as usize >= self.array.len() {
            return Ok((f32::NAN, (f32::NAN, f32::NAN)));
        }
        Ok((self.array[*x as usize][*y as usize], (0.0, 0.0)))
    }
}

#[allow(dead_code)]
impl ArrayDepth {
    pub(crate) fn new(array: Vec<Vec<f32>>) -> Self {
        ArrayDepth { array }
    }
}
