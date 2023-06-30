#[derive(Debug, thiserror::Error)]
pub(crate) enum Error {
   #[error("Undefined error")]
   /// The ode_solvers integrate method can return an error.
   Undefined,
   #[error("Argument passed was out of bounds")]
   /// The value k = |(kx, ky)| can only be positive.
   ArgumentOutOfBounds
}