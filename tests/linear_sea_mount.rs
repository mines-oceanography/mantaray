use mantaray;
use std::path::Path;

#[cfg(test)]
#[test]
/// verify the linear_sea_mount netcdf test works as expected. rays
/// should bend towards the circular island such that rays above the island
/// curve down and rays below the island curve up.
fn test_sea_mount() {
    let bathymetry_data = mantaray::bathymetry::cartesian::CartesianFile::new(Path::new(
        "tests/data/island_slice.nc",
    ));

    // top half
    let initial_waves: Vec<(f64, f64, f64, f64)> = (0..10)
        .map(|v| (-990.0, (v * 100) as f64, 0.01, 0.0))
        .collect();

    let waves = mantaray::ray::ManyRays::new(&bathymetry_data, None, &initial_waves);

    let results = waves.trace_many(0.0, 1000.0, 1.0);

    // verify the top half curves down towards the island. the x value
    // should be increasing. y value should stay the same or decrease.
    for res in results {
        match res {
            Some(res) => {
                let (_, data) = &res.get();
                // x values increase, y values decrease
                let mut last_x = data[0][0];
                let mut last_y = data[0][1];
                for r in data.iter() {
                    if r[0].is_nan() {
                        continue;
                    }
                    assert!(r[0] >= last_x);
                    assert!(r[1] <= last_y);
                    last_x = r[0];
                    last_y = r[1];
                }
            }
            None => (),
        }
    }

    // bottom half
    let initial_waves: Vec<(f64, f64, f64, f64)> = (1..10)
        .map(|v| (-990.0, (v * -100) as f64, 0.01, 0.0))
        .collect();

    let waves = mantaray::ray::ManyRays::new(&bathymetry_data, None, &initial_waves);

    let results = waves.trace_many(0.0, 1000.0, 1.0);

    // verify the bottom half curves up towards the island. x should
    // increase. y should increase.
    for res in &results {
        match res {
            Some(res) => {
                let (_, data) = &res.get();
                // x values increase, y values increase
                let mut last_x = data[0][0];
                let mut last_y = data[0][1];
                for r in data.iter() {
                    if r[0].is_nan() {
                        continue;
                    }
                    assert!(r[0] >= last_x);
                    assert!(r[1] >= last_y);
                    last_x = r[0];
                    last_y = r[1];
                }
            }
            None => (),
        }
    }
}
