//! Module containing interpolators
//! 
//! Contains the bilinear_interpolator function

#[allow(dead_code)]
/// Bilinear interpolation 
/// 
/// Performs operations to calculate bilinear interpolation at target point t
/// 
/// # Arguments
/// `points` : `&Vec<(i32, i32, f64)>`
/// - the known points with depth values. points must be in clockwise (relative)
///   order to each other with respect to the center of the square.
/// 
/// `target` : `&(f64, f64)`
/// - the target point must be contained within the square of the points.
/// 
/// # Panics
/// There are two ways for the function to panic, which should only happen due
/// to incorrect arguments passed.
/// - if the length of the inputted points vector is not 4.
/// - if the determinant is 0.
/// 
/// # Note
/// The points must be in correct order since the function assumes they are. It
/// will not give any error, but will return a value that is incorrect. In the
/// future, this function will enforce order of the points.
pub(crate) fn bilinear(points: &Vec<(i32, i32, f64)>, target: &(f64, f64)) -> f64 {
    // verify quadrilateral input
    assert!(points.len() == 4);

    // check if target is coincident with a point
    for point in points {
        if target.0 == point.0 as f64 && target.1 == point.1 as f64 {
            return point.2;
        }
    }

    // points are already in order
    let a = points[0];
    let b = points[1];
    let c = points[2];
    let d = points[3];

    // translate points and target with respect to a:
    let at = (0.0 , 0.0, a.2);
    let bt = (b.0 - a.0, b.1 - a.1, b.2);
    let ct = (c.0 - a.0, c.1 - a.1, c.2);
    let dt = (d.0 - a.0, d.1 - a.1, d.2);
    let tt =(target.0 - a.0 as f64, target.1 - a.1 as f64);
    println!("{:?}, {:?}, {:?}, {:?}", at, bt, ct, dt);

    // change basis of target point
    let det_bd = ((bt.0 * dt.1) - (dt.0 * bt.1)) as f64;
    assert!(det_bd != 0.0);
    // create inverse change of basis matrix
    let cbm = vec![
        vec![dt.1 as f64 / det_bd, -(dt.0 as f64 / det_bd)],
        vec![-(bt.1 as f64 / det_bd), bt.0 as f64 / det_bd]
    ];
    // calculate new target x and y coordinates (between 0 and 1)
    let x = cbm[0][0] * tt.0 + cbm[0][1] * tt.1;
    let y = cbm[1][0] * tt.0 + cbm[1][1] * tt.1;
    println!("x: {}, y: {}", x, y);

    // compute final value for the target position (bilinear interpolation)
    let a00 = a.2;
    let a10 = b.2 - a.2;  // change in the function's values at the points on the right and left at the same y
    let a01 = d.2 - a.2;  // change in the function's values at the points on the top and bottom at the same x
    let a11 = c.2 - a.2 - a10 - a01;  // change in x times the change in y

    a00 + a10 * x + a01 * y + a11 * x * y

}

#[test]
/// test single cases of the function against https://www.omnicalculator.com/math/bilinear-interpolation
fn test_interp() {
    // points must be in clockwise (relative) order to each other with respect to the center of the square.

    let check_interp = [
        (20.0, 23.0, -77, -19, 123, 145, 10.0, 20.0, 30.0, 40.0, 1230, 19.971951219512192),
    ];

    for (x, y, x1, y1, x2, y2, q11, q21, q12, q22, t, val) in check_interp {
        let points = vec![
            (x1 + t, y1 + t, q11),
            (x1 + t, y2 + t, q12),
            (x2 + t, y2 + t, q22),
            (x2 + t, y1 + t, q21),
        ];

        let target = (x + t as f64, y + t as f64);
        let ans = bilinear(&points, &target);
        assert!((ans - val).abs() < f64::EPSILON, "expected: {}. actual value: {}", val, ans);
    }
}

#[test]
/// test if the target is coincident with one of the input points
fn test_edges() {
    // x, y, x1, y1, x2, y2, q11, q21, q12, q22, t, val
    let check_interp = [
        (0.0, 0.0, 0, 0, 10, 10, 0.0, 5.0, 10.0, 15.0, 0, 0.0),
        (10.0, 0.0, 0, 0, 10, 10, 0.0, 5.0, 10.0, 15.0, 0, 5.0),
        (0.0, 10.0, 0, 0, 10, 10, 0.0, 5.0, 10.0, 15.0, 0, 10.0),
        (10.0, 10.0, 0, 0, 10, 10, 0.0, 5.0, 10.0, 15.0, 0, 15.0),
    ];

    for (x, y, x1, y1, x2, y2, q11, q21, q12, q22, t, val) in check_interp {
        let points = vec![
            (x1 + t, y1 + t, q11),
            (x1 + t, y2 + t, q12),
            (x2 + t, y2 + t, q22),
            (x2 + t, y1 + t, q21),
        ];

        let target = (x + t as f64, y + t as f64);
        let ans = bilinear(&points, &target);
        assert!((ans - val).abs() < f64::EPSILON, "expected: {}. actual value: {}", val, ans);
    }
}