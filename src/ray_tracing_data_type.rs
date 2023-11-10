//! Ray tracing data types and conversions
pub(crate) enum RayTracingDataType {
    VectorI32(Vec<i32>),
    VectorF32(Vec<f32>),
    VectorF64(Vec<f64>),
}

pub trait FromPrimitive {
    fn from_i32(n: i32) -> Option<Self>
    where
        Self: Sized;
    fn from_f32(n: f32) -> Option<Self>
    where
        Self: Sized;

    fn from_f64(n: f64) -> Option<Self>
    where
        Self: Sized;
}

impl FromPrimitive for i32 {
    fn from_i32(n: i32) -> Option<Self> {
        Some(n)
    }

    fn from_f32(n: f32) -> Option<Self>
    where
        Self: Sized,
    {
        Some(n as i32)
    }

    fn from_f64(n: f64) -> Option<Self>
    where
        Self: Sized,
    {
        Some(n as i32)
    }
}

impl FromPrimitive for f32 {
    fn from_i32(n: i32) -> Option<Self> {
        Some(n as f32)
    }

    fn from_f32(n: f32) -> Option<Self>
    where
        Self: Sized,
    {
        Some(n)
    }

    fn from_f64(n: f64) -> Option<Self>
    where
        Self: Sized,
    {
        Some(n as f32)
    }
}

impl FromPrimitive for f64 {
    fn from_i32(n: i32) -> Option<Self> {
        Some(n as f64)
    }

    fn from_f32(n: f32) -> Option<Self>
    where
        Self: Sized,
    {
        Some(n as f64)
    }

    fn from_f64(n: f64) -> Option<Self>
    where
        Self: Sized,
    {
        Some(n)
    }
}

pub fn convert_from_i32<T: FromPrimitive>(x: i32) -> Option<T> {
    T::from_i32(x)
}

pub fn convert_from_f32<T: FromPrimitive>(x: f32) -> Option<T> {
    T::from_f32(x)
}

pub fn convert_from_f64<T: FromPrimitive>(x: f64) -> Option<T> {
    T::from_f64(x)
}
