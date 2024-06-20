//! helper functions and types for integration testing

use mantaray::State;

pub const XINDEX: usize = 0;
pub const YINDEX: usize = 1;
pub const KX_INDEX: usize = 2;
pub const KY_INDEX: usize = 3;

/// assert that the value at the given index increases at each time step
pub fn assert_increase(data: &Vec<State>, index: usize) {
    let mut last = data[0][index];
    for r in data.iter().skip(1) {
        assert!(r[index] > last);
        last = r[index];
    }
}

/// assert that the value at the given index decreases at each time step
pub fn assert_decrease(data: &Vec<State>, index: usize) {
    let mut last = data[0][index];
    for r in data.iter().skip(1) {
        assert!(r[index] < last);
        last = r[index];
    }
}

/// assert that the value at the given index is the same at each time step
pub fn assert_same(data: &Vec<State>, index: usize) {
    let start_value = data[0][index];
    data.iter()
        .filter(|v| !v[0].is_nan())
        .for_each(|r| assert_eq!(r[index], start_value));
}
