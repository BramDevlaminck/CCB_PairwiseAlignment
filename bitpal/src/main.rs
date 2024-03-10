use clap::Parser;

use bitpal::bitpal;
use fasta_reader::read_fasta;

#[derive(Parser, Debug)]
#[clap(allow_negative_numbers = true)]
struct Args {
    /// The filename of the input fasta file containing the 2 sequences
    #[clap(short, long)]
    filename: String,
}

fn main() {
    let args = Args::parse();
    let Args { filename } = args;
    let (seq1, seq2) = read_fasta(&filename).expect("Failed to read the input fasta file");
    let seq1_chars = seq1.into_bytes();
    let seq2_chars = seq2.into_bytes();

    let result_score = bitpal(&seq1_chars, &seq2_chars);

    match result_score {
        Ok(score) => println!("The resulting match score is: {}", score),
        Err(error) => eprintln!("{}", error)
    }
}
