use std::error::Error;
use std::fmt::{Display, Formatter};

/// Error indicating that the provided input sequences to the BitPAl algorithm were too long
#[derive(Debug, PartialEq)]
pub struct InputTooLongError;

impl Display for InputTooLongError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "At least one of the provided sequences should be shorter than 64 characters")
    }
}

impl Error for InputTooLongError {}