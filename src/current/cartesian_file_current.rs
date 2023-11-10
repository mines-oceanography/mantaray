use std::path::Path;

use netcdf3::FileReader;

use crate::{
    error::Error,
    ray_tracing_data_type::{
        convert_from_f32, convert_from_f64, convert_from_i32, FromPrimitive, RayTracingDataType,
    },
};

const ENUM_VECTOR_I32: RayTracingDataType = RayTracingDataType::VectorI32(vec![]);
const ENUM_VECTOR_F32: RayTracingDataType = RayTracingDataType::VectorF32(vec![]);
const ENUM_VECTOR_F64: RayTracingDataType = RayTracingDataType::VectorF64(vec![]);

#[derive(Debug)]
/// Generic struct to accept different types of variables for input and save as
/// a consistent type.
pub(crate) struct CurrentCartesianFile<T> {
    vars_x_y_u_v: (Vec<T>, Vec<T>, Vec<T>, Vec<T>),
}

impl<T> CurrentCartesianFile<T>
where
    T: Sized + FromPrimitive,
{
    #[allow(dead_code)]

    pub(crate) fn new(
        path: &Path,
        x_type: &RayTracingDataType,
        y_type: &RayTracingDataType,
        u_type: &RayTracingDataType,
        v_type: &RayTracingDataType,
    ) -> Self {
        let mut data = FileReader::open(path).unwrap();

        let x_data = match x_type {
            RayTracingDataType::VectorI32(_) => data
                .read_var_i32("x")
                .unwrap()
                .iter()
                .map(|x| convert_from_i32(*x).unwrap())
                .collect(),
            RayTracingDataType::VectorF32(_) => data
                .read_var_f32("x")
                .unwrap()
                .iter()
                .map(|x| convert_from_f32(*x).unwrap())
                .collect(),
            RayTracingDataType::VectorF64(_) => data
                .read_var_f64("x")
                .unwrap()
                .iter()
                .map(|x| convert_from_f64(*x).unwrap())
                .collect(),
        };

        let y_data = match y_type {
            RayTracingDataType::VectorI32(_) => data
                .read_var_i32("y")
                .unwrap()
                .iter()
                .map(|x| convert_from_i32(*x).unwrap())
                .collect(),
            RayTracingDataType::VectorF32(_) => data
                .read_var_f32("y")
                .unwrap()
                .iter()
                .map(|x| convert_from_f32(*x).unwrap())
                .collect(),
            RayTracingDataType::VectorF64(_) => data
                .read_var_f64("y")
                .unwrap()
                .iter()
                .map(|x| convert_from_f64(*x).unwrap())
                .collect(),
        };

        let u_data = match u_type {
            RayTracingDataType::VectorI32(_) => data
                .read_var_i32("u")
                .unwrap()
                .iter()
                .map(|x| convert_from_i32(*x).unwrap())
                .collect(),
            RayTracingDataType::VectorF32(_) => data
                .read_var_f32("u")
                .unwrap()
                .iter()
                .map(|x| convert_from_f32(*x).unwrap())
                .collect(),
            RayTracingDataType::VectorF64(_) => data
                .read_var_f64("u")
                .unwrap()
                .iter()
                .map(|x| convert_from_f64(*x).unwrap())
                .collect(),
        };

        let v_data = match v_type {
            RayTracingDataType::VectorI32(_) => data
                .read_var_i32("v")
                .unwrap()
                .iter()
                .map(|x| convert_from_i32(*x).unwrap())
                .collect(),
            RayTracingDataType::VectorF32(_) => data
                .read_var_f32("v")
                .unwrap()
                .iter()
                .map(|x| convert_from_f32(*x).unwrap())
                .collect(),
            RayTracingDataType::VectorF64(_) => data
                .read_var_f64("v")
                .unwrap()
                .iter()
                .map(|x| convert_from_f64(*x).unwrap())
                .collect(),
        };

        let vars_x_y_u_v: (Vec<T>, Vec<T>, Vec<T>, Vec<T>) = (x_data, y_data, u_data, v_data);

        CurrentCartesianFile { vars_x_y_u_v }
    }
}

#[test]
fn test_generic() {
    let path = Path::new("dummy.nc");
    let current_data: CurrentCartesianFile<f32> = CurrentCartesianFile::new(
        &path,
        &ENUM_VECTOR_F64,
        &ENUM_VECTOR_F64,
        &ENUM_VECTOR_F64,
        &ENUM_VECTOR_F64,
    );
    dbg!(current_data);

    let current_data: CurrentCartesianFile<i32> = CurrentCartesianFile::new(
        &path,
        &ENUM_VECTOR_F64,
        &ENUM_VECTOR_F64,
        &ENUM_VECTOR_F64,
        &ENUM_VECTOR_F64,
    );
    dbg!(current_data);
}
