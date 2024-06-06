//! Trait with default implementations for converting an object that is
//! `Serializable` into a json string, writing itself, and saving itself in a
//! json file.

use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

use serde::Serialize;

use crate::error::Error;

/// Default implementations for converting an object that is `Serializable` into
/// a json string, writing itself, and saving itself in a json file.
pub trait WriteJson {
    /// Convert `Self` to a json string.
    ///
    /// # Returns
    ///
    /// json string of `Self`
    fn to_json_string(&self) -> String
    where
        Self: Serialize,
    {
        serde_json::to_string(&self).unwrap()
    }

    /// Write `Self` to a writer.
    ///
    /// # Arguments
    ///
    /// `writer` : `&mut W`
    /// - object that implements `Write`
    ///
    /// # Returns
    ///
    /// `Ok(usize)` : the number of bytes written
    ///
    /// `Err(Error)` : an error occurred while writing
    ///
    /// # Note
    ///
    /// This method writes `Self` as a json string.
    fn write_json<W: Write>(&self, writer: &mut W) -> Result<usize, Error>
    where
        Self: Serialize,
    {
        writer.write_all(self.to_json_string().as_bytes())?;
        writer.flush()?;
        Ok(self.to_json_string().as_bytes().len())
    }

    /// Save `Self` to a json file at the given path.
    ///
    /// # Arguments
    ///
    /// `path` : `&Path`
    /// - the path to save the file
    ///
    /// # Returns
    ///
    /// `Ok(usize)` : the number of bytes written
    ///
    /// `Err(Error)` : an error occurred while writing
    ///
    /// # Note
    ///
    /// This method writes `Self` as a json string at the given file path.
    fn save_json_file(&self, path: &Path) -> Result<usize, Error>
    where
        Self: Serialize,
    {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        self.write_json(&mut writer)
    }
}
