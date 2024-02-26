//! trait for writing to different file types
//!
//! Contains methods:
//! - write_to_csv_file

use std::path::Path;

use crate::error::Error;

pub(crate) trait Writable {
    /// write self to a csv file at given path location and with given delimiter
    fn write_to_csv_file(&self, path: &Path, delimiter: &str) -> Result<(), Error>;
}
