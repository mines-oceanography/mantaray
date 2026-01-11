//! Module containing interpolators
//!
//! Contains the bilinear_interpolator function

use tracing::trace;

use crate::error::{Error, Result};

#[derive(Debug)]
/// Handles a linear relationship
///
/// This was originally created to handle the dimensions of regular girds,
/// such latitude and longitude in a regularly spaced dataset, providing a
/// cheap conversion between the dimension, such as latitute, to the
/// correspondent index position.
pub(super) struct LinearFit<T> {
    /// Scale between the physical dimension and the index position.
    slope: T,
    /// Offset where the grid starts.
    intercept: T,
}

impl<T> LinearFit<T>
where
    T: Copy,
    T: std::ops::Sub<Output = T>,
    T: std::ops::Mul<Output = T>,
{
    #[allow(dead_code)]
    /// Convert to 'index' scale position
    ///
    /// Predict the index position of a given value. For instance, a result of
    /// `1.5` meants that the given value is between the second and third
    /// positions in the grid (first gridpoint is 0).
    pub(super) fn predict(&self, x: T) -> T {
        (x - self.intercept) * self.slope
    }
}

/*
#[cfg(test)]
mod test_linearfit {
    use super::*;

    #[test]
    fn test_fit_u64() {
        let lf = LinearFit::<u64> {
            slope: 2,
            intercept: 1,
        };
        assert_eq!(lf.predict(3), 7);
    }

    #[test]
    fn test_fit_f64() {
        let lf = LinearFit::<f64> {
            slope: 2.0,
            intercept: 1.0,
        };
        assert_eq!(lf.predict(3.0), 7.0);
    }
}
*/

/// Tolerance for linear relationship
///
/// If the ratio between the two dimensions deviate more than that, it will
/// not be considered a linear relationship.
const LINEAR_RELATION_TOLERANCE: f64 = 0.005;

impl LinearFit<f64> {
    /// Create a new LinearFit from a vector of values
    ///
    /// For a given vector of values, calculates the best linear relationship
    /// between the values and its index position with the purpose to estimate
    /// the index position of the closest value to a given target.
    ///
    /// This procedure also validates if a linear relationship is a good
    /// approximation with a threshold of 0.5% of tolerance.
    pub(super) fn from_fit(x: ndarray::ArrayD<f64>) -> Result<LinearFit<f64>> {
        let dx = &x.slice(ndarray::s![1..]) - &x.slice(ndarray::s![..-1]);
        let slope = dx.mean().expect("Failed to calculate mean");
        let criteria = ((dx - slope) / slope)
            .abs()
            .into_iter()
            .map(|v| v > LINEAR_RELATION_TOLERANCE)
            .any(|v| v);
        if criteria {
            return Err(Error::Undefined(
                "Linear is a bad approximation".to_string(),
            ));
        }
        Ok(LinearFit {
            slope,
            intercept: x[0],
        })
    }
}
