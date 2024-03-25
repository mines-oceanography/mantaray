//! # Data types

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
