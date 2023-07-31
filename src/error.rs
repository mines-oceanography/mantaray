//! Custom enum Error defining ray_tracing specific errors.
//! 
//! The enum depends on thiserror::Error.

#[derive(Debug, thiserror::Error)]
pub(crate) enum Error {
   #[error("Undefined error")]
   /// The undefined error will be used as a placeholder for all other errors.
   /// Two known errors that are currently undefined but need their own instance
   /// in the future are an index out of bounds error and an integration error.
   Undefined,
   #[error("Argument passed was out of bounds")]
   /// The value k = |(kx, ky)| can only be positive. If k <=0, the function will pass ArgumentOutOfBounds.
   ArgumentOutOfBounds,
   #[error("There was an error when computing the projection")]
   /// The projection should convert between latitude and longitude in degrees and x and y in meters.
   ProjectionError,
}