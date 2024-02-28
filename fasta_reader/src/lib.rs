use std::error::Error;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::path::{Path};

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn read_fasta(file: &str) -> Result<(String, String), Box<dyn Error>> {
    let mut sequences: Vec<String> = vec![];

    let mut current_sequence = String::new();
    for line in read_lines(file)?.flatten() {
        // skip header and push currently accumulated protein (if it is not empty, since that happens when, encountering the first header)
        if line.starts_with('>') {
            if !current_sequence.is_empty() {
                sequences.push(current_sequence);
                current_sequence = String::new();
            }
        } else {
            current_sequence += line.strip_suffix('\n').unwrap_or(&*line);
        }
    }
    sequences.push(current_sequence);

    let seq1 = sequences.swap_remove(0);
    let seq2 = sequences.swap_remove(0);

    Ok((seq1, seq2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_fasta_test() -> Result<(), Box<dyn Error>> {
        let (seq1, seq2) = read_fasta("../tests/input.fasta")?;
        assert_eq!(seq1, "GATTACA");
        assert_eq!(seq2, "GCATGCU");

        Ok(())
    }

    #[test]
    fn read_multiline_fasta_test() -> Result<(), Box<dyn Error>> {
        let (seq1, seq2) = read_fasta("../tests/multiline_input.fasta")?;
        assert_eq!(seq1, "GATTACA");
        assert_eq!(seq2, "GCATGCU");

        Ok(())
    }
}
