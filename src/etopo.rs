mod etopo {

    use std::collections::HashMap;
    
    fn main() {
    }

}

#[cfg(test)]
mod test_netcdf {

    #[test]
    fn test_main() {
        let etopo5 = netcdf::open("data.nc").expect("Could not open file");
    }
}