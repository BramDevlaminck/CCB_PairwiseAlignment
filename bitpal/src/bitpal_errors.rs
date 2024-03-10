use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Debug, PartialEq)]
pub struct InputTooLongError;

impl Display for InputTooLongError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "At least one of the provided sequences should be shorter than 64 characters")
    }
}

impl Error for InputTooLongError {}