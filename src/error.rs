//! Custom enum Error defining ray_tracing specific errors.
//!
//! The enum depends on thiserror::Error.

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("Argument passed was out of bounds")]
    /// The value k = |(kx, ky)| can only be positive. If k <=0, the function will pass ArgumentOutOfBounds.
    ArgumentOutOfBounds,

    #[error("One or more of the points surrounding the nearest point are out of bounds.")]
    /// This error is returned when `nearest_point` returns a point in bounds,
    /// but `four_corners` is still out of bounds.
    CornersOutOfBounds,

    #[error("Argument passed was not a valid option")]
    /// The argument passed was not a valid option
    InvalidArgument,

    #[error("Index passed was out of bounds")]
    /// The index is out of bounds of the array and would panic if attempted to
    /// access array.
    IndexOutOfBounds,

    #[allow(clippy::enum_variant_names)] // tell clippy the name is ok
    #[error(transparent)]
    // IO error from std::io
    IOError(#[from] std::io::Error),

    #[allow(clippy::enum_variant_names)] // tell clippy the name is ok
    #[error(transparent)]
    // Integration error from ode_solvers
    IntegrationError(#[from] ode_solvers::dop_shared::IntegrationError),

    #[error("The requested point has no nearest point")]
    /// The target point was either outside the domain or closest to the edge of
    /// the domain.
    NoNearestPoint,
}

pub type Result<T> = core::result::Result<T, Error>;
