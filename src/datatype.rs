//! # Data types

use crate::wave_ray_path::State;

#[derive(Clone, Debug)]
/// A point in 2D cartesian space
///
/// A `Point` is composed by `x` and `y`, expected to be in meters.
pub(crate) struct Point<T> {
    x: T,
    y: T,
}

#[allow(dead_code)]
impl<T> Point<T> {
    /// Create a new `Point` with the given `x` and `y` coordinates.
    ///
    pub(crate) fn new(x: T, y: T) -> Self {
        Point { x, y }
    }

    /// Get the `x` coordinate of the `Point`.
    ///
    pub(crate) fn x(&self) -> &T {
        &self.x
    }

    /// Get the `y` coordinate of the `Point`.
    ///
    pub(crate) fn y(&self) -> &T {
        &self.y
    }
}

/// A 2D geolocation in a 2D space
///
/// A `Coordinate` is composed by `lat` and `lon`, expected to be in decimal
/// degrees. For instance, the latitude of the North Pole is 90, and a
/// latitude of -10.5 is equivalent to 10 degrees and 30 minutes South.
pub(crate) struct Coordinate<T> {
    lat: T,
    lon: T,
}

#[allow(dead_code)]
impl<T> Coordinate<T> {
    /// Create a new `Coordinate` with the given `lat` and `lon` coordinates.
    ///
    fn new(lat: T, lon: T) -> Self {
        Coordinate { lat, lon }
    }

    fn lat(&self) -> &T {
        &self.lat
    }

    fn lon(&self) -> &T {
        &self.lon
    }
}

#[allow(dead_code)]
#[derive(Clone, Debug, PartialEq)]
/// The current in a 2D cartesian point
///
/// A `Current` is composed by `u` and `v`, expected to be in meters per
/// second.
pub(crate) struct Current<T> {
    u: T,
    v: T,
}

#[allow(dead_code)]
impl<T> Current<T> {
    pub(crate) fn new(u: T, v: T) -> Self {
        Current { u, v }
    }

    pub(crate) fn u(&self) -> &T {
        &self.u
    }

    pub(crate) fn v(&self) -> &T {
        &self.v
    }
}

#[allow(dead_code)]
#[derive(Clone, Debug)]
/// A wave number in 2D cartesian space
pub(crate) struct WaveNumber<T> {
    kx: T,
    ky: T,
}

#[allow(dead_code)]
impl<T> WaveNumber<T> {
    /// create a new wave number from the given `kx` and `ky` values
    pub(crate) fn new(kx: T, ky: T) -> Self {
        WaveNumber { kx, ky }
    }

    /// get the x component of the wave number
    pub(crate) fn kx(&self) -> &T {
        &self.kx
    }

    /// get the y component of the wave number
    pub(crate) fn ky(&self) -> &T {
        &self.ky
    }
}

#[allow(dead_code)]
#[derive(Clone, Debug)]
/// a ray state is the point and wave number of the ray
pub(crate) struct RayState<T> {
    // Position in 2D cartesian space.
    point: Point<T>,
    // Wave number in 2D cartesian space.
    wave_number: WaveNumber<T>,
}

impl<T> RayState<T> {
    /// create a new `RayState`
    pub(crate) fn new(point: Point<T>, wave_number: WaveNumber<T>) -> Self {
        RayState { point, wave_number }
    }

    fn point(&self) -> &Point<T> {
        &self.point
    }

    /// get the wave number of the ray state
    pub(crate) fn wave_number(&self) -> &WaveNumber<T> {
        &self.wave_number
    }
}

impl From<RayState<f64>> for State {
    /// convert mantaray's `RayState` into `State` object used by `ode_solvers`.
    fn from(value: RayState<f64>) -> Self {
        State::new(
            *value.point().x(),
            *value.point().y(),
            *value.wave_number().kx(),
            *value.wave_number().ky(),
        )
    }
}

// Possible names:
// - RayPath
// - RayTrajectory
// - Beam
//

#[derive(Clone, Debug)]
/// A single wave ray in 2D cartesian space
///
/// Parameters defining the evolution of an wave ray.
///
/// Note that we might generalize later to allow coordinate as an alternative
/// to point.
pub(crate) struct Ray<T> {
    // Relative time in seconds. Initial condition is t=0.
    time: Vec<f32>,
    state: Vec<RayState<T>>,
    // Depth in meters.
    depth: Vec<f32>,
    // Current in 2D cartesian space.
    current: Vec<Current<T>>,
}

#[allow(dead_code)]
impl<T> Ray<T> {
    /// Create a new `Ray`
    fn new() -> Self {
        Ray {
            time: Vec::new(),
            state: Vec::new(),
            depth: Vec::new(),
            current: Vec::new(),
        }
    }

    /// Push a new state and environment information to the `Ray`
    fn push(&mut self, time: f32, state: RayState<T>, depth: f32, current: Current<T>) {
        self.time.push(time);
        self.state.push(state);
        self.depth.push(depth);
        self.current.push(current);
    }
}

#[cfg(test)]
mod test_ray {
    use super::*;

    #[test]
    fn test_new() {
        let ray: Ray<f32> = Ray::new();
        assert_eq!(ray.time.len(), 0);
        assert_eq!(ray.state.len(), 0);
        assert_eq!(ray.depth.len(), 0);
        assert_eq!(ray.current.len(), 0);
    }

    #[test]
    fn test_push() {
        let mut ray: Ray<f32> = Ray::new();
        let point = Point::new(1.0, 2.0);
        let wave_number = WaveNumber { kx: 3.0, ky: 4.0 };
        let current = Current::new(5.0, 6.0);
        ray.push(0.0, RayState { point, wave_number }, 100.0, current);
        assert_eq!(ray.time.len(), 1);
        assert_eq!(ray.state.len(), 1);
        assert_eq!(ray.depth.len(), 1);
        assert_eq!(ray.current.len(), 1);
    }
}

pub(crate) struct Bundle<T> {
    rays: Vec<Ray<T>>,
}

#[allow(dead_code)]
impl<T> Bundle<T> {
    fn new() -> Self {
        Bundle { rays: Vec::new() }
    }

    fn push(&mut self, ray: Ray<T>) {
        self.rays.push(ray);
    }
}

#[derive(Debug, PartialEq)]
/// A gradient in 2D space
pub(crate) struct Gradient<T> {
    dx: T,
    dy: T,
}

impl<T> Gradient<T> {
    pub(crate) fn new(dx: T, dy: T) -> Self {
        Gradient { dx, dy }
    }

    pub(crate) fn dx(&self) -> &T {
        &self.dx
    }

    pub(crate) fn dy(&self) -> &T {
        &self.dy
    }
}
