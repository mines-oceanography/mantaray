//! Module containing interpolators
//!
//! Contains the bilinear_interpolator function

use std::collections::HashMap;

use tracing::trace;

use crate::datatype::Point;
use crate::error::{Error, Result};
use crate::io::Dataset;

#[allow(dead_code)]
/// Bilinear interpolation
///
/// Performs operations to calculate bilinear interpolation at target point
///
/// # Arguments
/// `points` : `&Vec<(i32, i32, f64)>`
/// - the known points with depth values. points must be in clockwise (relative)
///   order to each other with respect to the center of the square.
///
/// `target` : `&(f64, f64)`
/// - the target point must be contained within the square of the points.
///
/// # Returns
/// `Result<f32, Error>`
/// - `Ok(f32)` : interpolated depth
/// - `Err(Error)` : argument passed `points` is invalid
///
/// # Errors
/// `Error::InvalidArgument` : either the number of points is not equal to 4, or
/// the determinant of the change of basis matrix equals zero.
///
/// # Note
/// The points must be in correct order since the function assumes they are. It
/// will not give any error, but will return a value that is incorrect. In the
/// future, this function will enforce order of the points.
pub(crate) fn bilinear(points: &Vec<(f32, f32, f32)>, target: &(f32, f32)) -> Result<f32> {
    // verify quadrilateral input
    if points.len() != 4 {
        return Err(Error::InvalidArgument);
    }

    // check if target is coincident with a point
    for point in points {
        if target.0 == point.0 && target.1 == point.1 {
            return Ok(point.2);
        }
    }

    // points are already in order
    let a = points[0];
    let b = points[1];
    let c = points[2];
    let d = points[3];

    // translate points b, d, and target with respect to a:
    let bt = (b.0 - a.0, b.1 - a.1, b.2);
    let dt = (d.0 - a.0, d.1 - a.1, d.2);
    let tt = (target.0 - a.0, target.1 - a.1);

    // change basis of target point
    let det_bd = (bt.0 * dt.1) - (dt.0 * bt.1);
    if det_bd == 0.0 {
        return Err(Error::InvalidArgument);
    }
    // create inverse change of basis matrix
    let cbm = vec![
        vec![dt.1 / det_bd, -(dt.0 / det_bd)],
        vec![-(bt.1 / det_bd), bt.0 / det_bd],
    ];
    // calculate new target x and y coordinates (between 0 and 1)
    let x = cbm[0][0] * tt.0 + cbm[0][1] * tt.1;
    let y = cbm[1][0] * tt.0 + cbm[1][1] * tt.1;

    // compute final value for the target position (bilinear interpolation)
    let a00 = a.2;
    let a10 = b.2 - a.2; // change in the function's values at the points on the right and left at the same y
    let a01 = d.2 - a.2; // change in the function's values at the points on the top and bottom at the same x
    let a11 = c.2 - a.2 - a10 - a01; // change in x times the change in y

    Ok(a00 + a10 * x + a01 * y + a11 * x * y)
}

#[test]
/// test single cases of the function against https://www.omnicalculator.com/math/bilinear-interpolation
fn test_interp() {
    // points must be in clockwise (relative) order to each other with respect to the center of the square.

    let check_interp = [(
        20.0,
        23.0,
        -77.0,
        -19.0,
        123.0,
        145.0,
        10.0,
        20.0,
        30.0,
        40.0,
        1230.0,
        19.971951219512192,
    )];

    for (x, y, x1, y1, x2, y2, q11, q21, q12, q22, t, val) in check_interp {
        let points = vec![
            (x1 + t, y1 + t, q11),
            (x1 + t, y2 + t, q12),
            (x2 + t, y2 + t, q22),
            (x2 + t, y1 + t, q21),
        ];

        let target = (x + t, y + t);
        let ans = bilinear(&points, &target).unwrap();
        assert!(
            (ans - val).abs() < f32::EPSILON,
            "expected: {}. actual value: {}",
            val,
            ans
        );
    }
}

#[test]
/// test if the target is coincident with one of the input points
fn test_edges() {
    // x, y, x1, y1, x2, y2, q11, q21, q12, q22, t, val
    let check_interp = [
        (
            0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 0.0, 5.0, 10.0, 15.0, 0.0, 0.0,
        ),
        (
            10.0, 0.0, 0.0, 0.0, 10.0, 10.0, 0.0, 5.0, 10.0, 15.0, 0.0, 5.0,
        ),
        (
            0.0, 10.0, 0.0, 0.0, 10.0, 10.0, 0.0, 5.0, 10.0, 15.0, 0.0, 10.0,
        ),
        (
            10.0, 10.0, 0.0, 0.0, 10.0, 10.0, 0.0, 5.0, 10.0, 15.0, 0.0, 15.0,
        ),
    ];

    for (x, y, x1, y1, x2, y2, q11, q21, q12, q22, t, val) in check_interp {
        let points = vec![
            (x1 + t, y1 + t, q11),
            (x1 + t, y2 + t, q12),
            (x2 + t, y2 + t, q22),
            (x2 + t, y1 + t, q21),
        ];

        let target = (x + t, y + t);
        let ans = bilinear(&points, &target).unwrap();
        assert!(
            (ans - val).abs() < f32::EPSILON,
            "expected: {}. actual value: {}",
            val,
            ans
        );
    }
}

#[derive(Debug)]
/// Handles a linear relationship
///
/// This was originally created to handle the dimensions of regular girds,
/// such latitude and longitude in a regularly spaced dataset, providing a
/// cheap conversion between the dimension, such as latitute, to the
/// correspondent index position.
struct LinearFit<T> {
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
    fn predict(&self, x: T) -> T {
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
    fn from_fit(x: ndarray::ArrayD<f64>) -> Result<LinearFit<f64>> {
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

struct RegularGrid<'a> {
    // dataset: &'a dyn Dataset,
    dataset: Box<dyn Dataset + 'a>,
    x_size: usize,
    x_map: LinearFit<f64>,
    y_size: usize,
    y_map: LinearFit<f64>,
    // Save the dimensions order: ij or ji
    dimension_order: HashMap<String, String>,
}

impl<'a> RegularGrid<'a> {
    /*
        fn validate_as_regular_grid() {
            // Open a full 1D array for x and another for y
            // calculate delta and confirm that all are <0.01 diference (criteria for linear)
            //
        }
    */
    #[allow(dead_code)]
    fn open(dataset: impl Dataset + 'a, varname_x: &str, varname_y: &str) -> Result<Self> {
        // confirm it is linear
        // Define A & B coefficients
        // get i_size and j_size

        // let dataset = netcdf::open(&file).unwrap();

        // Identify the variables that have the user defined dimensions
        // and create a map on the dimenson order

        let dimension_order = dataset.dimensions_order(varname_x, varname_y);

        let x_size = dataset.dimension_len(varname_x).unwrap();

        let x_map = LinearFit::from_fit(dataset.values(varname_x).unwrap())?;

        let y_size = dataset.dimension_len(varname_y).unwrap();
        let y_map = LinearFit::from_fit(dataset.values(varname_x).unwrap())?;

        Ok(Self {
            // dataset: dataset,
            dataset: Box::new(dataset),
            x_size,
            x_map,
            y_size,
            y_map,
            dimension_order,
        })
    }

    #[allow(dead_code)]
    /// Get the nearest `varname` value to the given `x` and `y` coordinates
    fn nearest(&self, varname: &str, point: Point<f64>) -> Result<f32> {
        // Error message if less than 0
        let i = self.x_map.predict(*point.x()).round() as usize;
        if i >= self.x_size {
            return Err(Error::IndexOutOfBounds);
        }
        let j = self.y_map.predict(*point.y()).round() as usize;
        if j >= self.y_size {
            return Err(Error::IndexOutOfBounds);
        }

        match self.dimension_order.get(varname) {
            Some(v) => match v.as_str() {
                "xy" => {
                    trace!("Assuming dimension order is 'xy'");
                    self.dataset.get_variable(varname, i, j)
                }
                "yx" => {
                    trace!("Assuming dimension order is 'yx'");
                    self.dataset.get_variable(varname, j, i)
                }
                _ => return Err(Error::Undefined("Dimension order not found".to_string())),
            },
            _ => return Err(Error::Undefined("Variable not found".to_string())),
        }
    }
}

/*

impl<T: num_traits::float::Float> RegularGrid<T> {
impl<T: num_traits::float::FloatCore> RegularGrid<T> {
impl<T: num_traits::real::Real> RegularGrid<T> {


impl<T> RegularGrid<T> {
    fn nearest(&self, x: T, y: T) -> Result<(usize, usize)>
    where
        T: Clone,
        T: std::ops::Mul<Output = T>,
        T: std::ops::Add<Output = T>,
        f64: From<T>,
        // usize: From<T>,
    {
        let x = self.x.fit(x);
        let y = self.y.fit(y);

        /*
        if x < 0).unwrap() || x >= T::from(self.x_size).unwrap() {
            return Err(Error::InvalidArgument);
        }

        if y < T::from(0).unwrap() || y >= T::from(self.y_size).unwrap() {
            return Err(Error::InvalidArgument);
        }
        */

        let x = f64::round(x.into());
        let y = f64::round(y.into());

        Ok((x as usize, y as usize))
    }
    /*
    impl<T: std::ops::Mul> RegularGrid<T> {
        fn nearest(&self, x: T, y: T) -> Result<(usize, usize)>
        where
            <T as std::ops::Mul>::Output: std::ops::Add<T>,
            usize: From<<<T as std::ops::Mul>::Output as std::ops::Add<T>>::Output>,
        {
            let x = self.x.slope * x + self.x.intercept;
            let y = self.y.slope * y + self.y.intercept;

            let x: usize = x.try_into().unwrap();
            let y: usize = y.try_into().unwrap();

            /*
            if x < 0 | x >= self.x_size as T {
                return Err(Error::InvalidArgument);
            }
            if y < 0 | y >= self.y_size as T {
                return Err(Error::InvalidArgument);
            }

            let x = x.round();
            let y = y.round();

            if x < T::from(0).unwrap() || x >= T::from(self.x_size).unwrap() {
                return Err(Error::InvalidArgument);
            }

            if y < T::from(0).unwrap() || y >= T::from(self.y_size).unwrap() {
                return Err(Error::InvalidArgument);
            }
            */

            Ok((x, y))
        }
    }
    */
}

#[cfg(test)]
mod test_regulargrid {
    use super::*;

    #[test]
    fn test_nearest() {
        let rg = RegularGrid {
            x: LinearFit {
                slope: 1.0,
                intercept: 0.0,
            },
            x_size: 10,
            y: LinearFit {
                slope: 1.0,
                intercept: 0.0,
            },
            y_size: 10,
        };

        assert_eq!(rg.nearest(3.0, 3.0).unwrap(), (3, 3));
        assert_eq!(rg.nearest(3.1, 3.1).unwrap(), (3, 3));
        assert_eq!(rg.nearest(3.5, 3.5).unwrap(), (4, 4));
        assert_eq!(rg.nearest(3.9, 3.9).unwrap(), (4, 4));
    }
}
*/
