//! Custom enum Error defining ray_tracing specific errors.
//!
//! The enum depends on thiserror::Error.

#[derive(Debug, thiserror::Error)]
#[allow(clippy::enum_variant_names)] // tell clippy the name is ok
pub(crate) enum Error {
    #[error("Argument passed was out of bounds")]
    /// The value k = |(kx, ky)| can only be positive. If k <=0, the function will pass ArgumentOutOfBounds.
    ArgumentOutOfBounds,

    #[error("Argument passed was not a valid option")]
    /// The argument passed was not a valid option
    InvalidArgument,

    #[error("Index passed was out of bounds")]
    /// The index is out of bounds of the array and would panic if attempted to
    /// access array.
    IndexOutOfBounds,

    #[error("Generic error: {0}")]
    /// Temporary error type. Any undefined error should be eventually
    /// replaced by a permanent type.
    Undefined(String),

    #[error(transparent)]
    // IO error from std::io
    IOError(#[from] std::io::Error),

    #[error(transparent)]
    // Integration error from ode_solvers
    IntegrationError(#[from] ode_solvers::dop_shared::IntegrationError),

    #[error(transparent)]
    // ReadError from netcdf3
    ReadError(#[from] netcdf3::error::ReadError),
}

pub(crate) type Result<T> = core::result::Result<T, Error>;
