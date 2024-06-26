//! helper functions and types for integration testing

use mantaray::State;

pub const XINDEX: usize = 0;
pub const YINDEX: usize = 1;
pub const KX_INDEX: usize = 2;
pub const KY_INDEX: usize = 3;

/// true if the value at the given index increases at each time step
pub fn increase(data: &Vec<State>, index: usize) -> bool {
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
pub fn decrease(data: &Vec<State>, index: usize) -> bool {
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
pub fn same(data: &Vec<State>, index: usize) -> bool {
    let mut last = data[0][index];
    for r in data.iter().filter(|v| !v[0].is_nan()).skip(1) {
        if !(r[index] == last) {
            return false;
        }
        last = r[index];
    }
    return true;
}

/// assert that the value at the given index stays the same until the condition.
/// After that point, the value can be greater than or equal to the previous
/// value.
pub fn assert_can_increase_after(data: &Vec<State>, index: usize, condition: fn(&State) -> bool) {
    let mut last_value = data[0][index];
    for state in data.iter().filter(|v| !v[0].is_nan()).skip(1) {
        let value = state[KX_INDEX];
        if !condition(state) {
            assert_eq!(last_value, value);
            last_value = value;
        } else {
            assert!(value >= last_value);
        }
    }
}
