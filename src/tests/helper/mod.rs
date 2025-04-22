//! helper functions and types for integration testing

use crate::State;

pub(crate) const XINDEX: usize = 0;
pub(crate) const YINDEX: usize = 1;
pub(crate) const KX_INDEX: usize = 2;
pub(crate) const KY_INDEX: usize = 3;

#[cfg(test)]

/// true if the value at the given index increases at each time step
pub(crate) fn increase(data: &Vec<State>, index: usize) -> bool {
    let mut last = data[0][index];
    for r in data.iter().filter(|v| !v[0].is_nan()).skip(1) {
        if !(r[index] > last) {
            return false;
        }
        last = r[index];
    }
    return true;
}

/// true if the value at the given index decreases at each time step
pub(crate) fn decrease(data: &Vec<State>, index: usize) -> bool {
    let mut last = data[0][index];
    for r in data.iter().filter(|v| !v[0].is_nan()).skip(1) {
        if !(r[index] < last) {
            return false;
        }
        last = r[index];
    }
    return true;
}

/// true if the value at the given index is exactly the same at each time step
pub(crate) fn same(data: &Vec<State>, index: usize) -> bool {
    let mut last = data[0][index];
    for r in data.iter().filter(|v| !v[0].is_nan()).skip(1) {
        if !(r[index] == last) {
            println!("Expected {last} but got {}", r[index]);
            return false;
        }
        last = r[index];
    }
    return true;
}
