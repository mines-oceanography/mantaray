//! # Data types

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
    fn new(x: T, y: T) -> Self {
        Point { x, y }
    }

    /// Get the `x` coordinate of the `Point`.
    ///
    fn x(&self) -> &T {
        &self.x
    }

    /// Get the `y` coordinate of the `Point`.
    ///
    fn y(&self) -> &T {
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
#[derive(Clone, Debug)]
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
    fn new(u: T, v: T) -> Self {
        Current { u, v }
    }

    fn u(&self) -> &T {
        &self.u
    }

    fn v(&self) -> &T {
        &self.v
    }
}

#[derive(Clone, Debug)]
/// A wave number in 2D cartesian space
struct WaveNumber<T> {
    kx: T,
    ky: T,
}

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
    // Position in 2D cartesian space.
    point: Vec<Point<T>>,
    // Wave number in 2D cartesian space.
    wave_number: Vec<WaveNumber<T>>,
    // Depth in meters.
    depth: Vec<f32>,
    // Current in 2D cartesian space.
    current: Vec<Current<T>>,
}

#[allow(dead_code)]
impl<T> Ray<T> {
    fn new() -> Self {
        Ray {
            time: Vec::new(),
            point: Vec::new(),
            wave_number: Vec::new(),
            depth: Vec::new(),
            current: Vec::new(),
        }
    }

    fn push(
        &mut self,
        time: f32,
        point: Point<T>,
        wave_number: WaveNumber<T>,
        depth: f32,
        current: Current<T>,
    ) {
        self.time.push(time);
        self.point.push(point);
        self.wave_number.push(wave_number);
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
        assert_eq!(ray.point.len(), 0);
        assert_eq!(ray.wave_number.len(), 0);
        assert_eq!(ray.depth.len(), 0);
        assert_eq!(ray.current.len(), 0);
    }

    #[test]
    fn test_push() {
        let mut ray: Ray<f32> = Ray::new();
        let point = Point::new(1.0, 2.0);
        let wave_number = WaveNumber { kx: 3.0, ky: 4.0 };
        let current = Current::new(5.0, 6.0);
        ray.push(0.0, point, wave_number, 7.0, current);
        assert_eq!(ray.time.len(), 1);
        assert_eq!(ray.point.len(), 1);
        assert_eq!(ray.wave_number.len(), 1);
        assert_eq!(ray.depth.len(), 1);
        assert_eq!(ray.current.len(), 1);
    }
}

pub(crate) struct RayBundle<T> {
    rays: Vec<Ray<T>>,
}

impl<T> RayBundle<T> {
    fn new() -> Self {
        RayBundle { rays: Vec::new() }
    }

    fn push(&mut self, ray: Ray<T>) {
        self.rays.push(ray);
    }
}
