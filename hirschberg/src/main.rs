use std::error::Error;
use std::str::from_utf8;
use clap::Parser;
use fasta_reader::read_fasta;
use hirschberg::hirschberg;

/// parse the arguments
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

fn main() -> Result<(), Box<dyn Error>>{
    let args = Args::parse();
    let Args { filename, match_score, mismatch_score, gap_score } = args;
    let (seq1, seq2) = read_fasta(&filename)?;
    let seq1_chars = seq1.into_bytes();
    let seq2_chars = seq2.into_bytes();

    let (aligned_seq1, diff_line, aligned_seq2) = hirschberg(&seq1_chars, &seq2_chars, match_score, mismatch_score, gap_score)?;

    println!("Aligned sequences:");
    println!("{}", from_utf8(&*aligned_seq1)?);
    println!("{}", from_utf8(&*diff_line)?);
    println!("{}", from_utf8(&*aligned_seq2)?);

    Ok(())
}
