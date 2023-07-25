//! BathymetryData trait

/// A trait used to give the function get_depth
pub(crate) trait BathymetryData {
    /// Returns the nearest depth for the given lat, lon coordinate.
    fn get_depth(&self, x: &f32, y: &f32) -> f64;
}
