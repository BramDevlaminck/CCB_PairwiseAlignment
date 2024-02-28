use std::error::Error;
use clap::Parser;
use fasta_reader::read_fasta;
use needleman_wunsch_delta_encoded::construct_delta_matrices;

#[derive(Parser, Debug)]
#[clap(allow_negative_numbers = true)]
struct Args {
    /// The input file name
    #[clap(short, long)]
    filename: String,
    /// The score used when there is a match
    #[clap(short, long, default_value_t = 1)]
    match_score: i32,
    /// The score used when there is a mismatch
    #[clap(short = 'i', long, default_value_t = - 1)]
    mismatch_score: i32,
    /// The score used when there is a gap
    #[clap(short, long, default_value_t = - 3)]
    gap_score: i32,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let Args { filename, match_score, mismatch_score, gap_score } = args;
    let (seq1, seq2) = read_fasta(&filename)?;
    let seq1_chars = seq1.into_bytes();
    let seq2_chars = seq2.into_bytes();

    let (delta_h, delta_v) = construct_delta_matrices(&seq1_chars, &seq2_chars, match_score, mismatch_score, gap_score);

    let row = seq2_chars.len();
    let mut delta_h_score = row as i32 * gap_score;
    for col in 1..=seq1_chars.len() {
        delta_h_score += delta_h[row][col].unwrap();
    }

    let col = seq1_chars.len();
    let mut delta_v_score = col as i32 * gap_score;
    for row in 1..=seq2_chars.len() {
        delta_v_score += delta_v[row][col].unwrap();
    }

    println!("score according to the delta h matrix {}", delta_h_score);
    println!("score according to the delta v matrix {}", delta_v_score);

    Ok(())
}