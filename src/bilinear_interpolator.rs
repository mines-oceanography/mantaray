//! Bilinear interpolator
//! 
//! Contains the bilinear_interpolator function

/// Bilinear interpolation 
/// 
/// Performs operations to calculate bilinear interpolation at target point t
/// 
/// # Arguments
/// `points` : `&mut Vec<(i32, i32, f64)>`
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
fn bilinear_interpolator(points: &mut Vec<(i32, i32, f64)>, target: &(f64, f64)) -> f64 {
    // verify quadrilateral input
    assert!(points.len() == 4);

    // points are already in order
    let a = points[0];
    let b = points[1];
    let c = points[2];
    let d = points[3];

    // commented below: to order the points from a random set does NOT work! maybe check logic again later:
    // let a = points.remove(0);
    // let b = *points.iter().min_by_key(|p| (p.0 - a.0).pow(2) + (p.1 - a.1).pow(2)).unwrap();
    // let c = *points.iter().max_by_key(|p| (p.0 - a.0).pow(2) + (p.1 - a.1).pow(2)).unwrap();
    // let d = *points.iter().find(|p| !(p.0 == a.0 && p.1 == a.1) && !(p.0 == b.0 && p.1 == b.1) && !(p.0 == c.0 && p.1 == c.1)).unwrap();
    // println!("{:?}, {:?}, {:?}, {:?}", a, b, c, d);

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
    let q11 = 10.0;
    let q21 = -10.0;
    let q12 = -10.0;
    let q22 = 10.0;

    let mut points = vec![
        (0, 0, q11),
        (5, 5, q21),
        (10, 0, q22),
        (5, -5, q12),
    ];
    let target = (5.0, 0.0);
    let ans = bilinear_interpolator(&mut  points, &target);
    assert!((ans - 0.0).abs() < f64::EPSILON, "actual value: {}", ans);
}