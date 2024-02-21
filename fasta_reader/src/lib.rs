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

    for (i, line) in read_lines(file)?.flatten().enumerate() {
        // even lines only contain headers, skip them
        if i % 2 == 0 {
            continue;
        }
        sequences.push(line);
    }

    let seq1 = sequences.swap_remove(0);
    let seq2 = sequences.swap_remove(0);

    Ok((seq1, seq2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_fasta_test() -> Result<(), Box<dyn Error>> {
        let (seq1, seq2) = read_fasta( "../tests/input.fasta")?;
        assert_eq!(seq1, "GATTACA");
        assert_eq!(seq2, "GCATGCU");

        Ok(())
    }
}
