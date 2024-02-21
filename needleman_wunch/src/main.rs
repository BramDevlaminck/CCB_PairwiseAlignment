use std::cmp::max;
use std::error::Error;
use std::str::from_utf8;
use clap::Parser;

/// parse the arguments
#[derive(Parser, Debug)]
#[clap(allow_negative_numbers = true)]
struct Args {
    /// Sequence 1
    seq1: String,
    /// Sequence 2
    seq2: String,
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
    let Args { seq1, seq2, match_score, mismatch_score, gap_score } = args;
    let seq1_chars = seq1.into_bytes();
    let seq2_chars = seq2.into_bytes();

    let mut matrix: Vec<Vec<i32>> = vec![vec![0; seq1_chars.len() + 1]; seq2_chars.len() + 1];

    // initialize the first column and first row
    matrix[0] = (0..seq1_chars.len() + 1).map(|i| i as i32 * gap_score).collect();
    for (index, init_value) in (0..seq2_chars.len() + 1).enumerate() {
        matrix[index][0] = init_value as i32 * gap_score;
    }

    // fill in the matrix
    for row in 1..seq2_chars.len() + 1 {
        let current_seq2_char = seq2_chars[row - 1];
        for col in 1..seq1_chars.len() + 1 {
            let diag_score = matrix[row - 1][col - 1] + if seq1_chars[col - 1] == current_seq2_char { match_score } else { mismatch_score };
            matrix[row][col] = max(
                max(diag_score, matrix[row][col - 1] + gap_score),
                matrix[row - 1][col] + gap_score
            );
        }
    }

    let mut current_row = seq2_chars.len();
    let mut current_col = seq1_chars.len();
    println!("Global score is: {}", matrix[current_row][current_col]);

    // backtrack to find alignment
    let mut aligned_seq1: Vec<u8> = vec![];
    let mut aligned_seq2: Vec<u8> = vec![];
    let mut diff_line: Vec<u8> = vec![];

    while current_col != 0 && current_row != 0 {
        if current_col != 0 && matrix[current_row][current_col] == matrix[current_row][current_col - 1] + gap_score {
            aligned_seq1.insert(0, seq1_chars[current_col - 1]);
            diff_line.insert(0, b' ');
            aligned_seq2.insert(0, b'-');
            current_col -= 1;
            continue
        }

        if current_row != 0 && matrix[current_row][current_col] == matrix[current_row - 1][current_col] + gap_score {
            aligned_seq1.insert(0, b'-');
            diff_line.insert(0, b' ');
            aligned_seq2.insert(0, seq2_chars[current_row - 1]);
            current_row -= 1;
            continue
        }

        // diagonal case
        aligned_seq1.insert(0, seq1_chars[current_col - 1]);
        aligned_seq2.insert(0, seq2_chars[current_row - 1]);
        let char_to_insert = if seq1_chars[current_col - 1] == seq2_chars[current_row - 1] { b'|' } else { b'*' };
        diff_line.insert(0, char_to_insert);
        current_row -= 1;
        current_col -= 1;
    }

    println!();
    println!("Aligned sequences:");
    println!("{}", from_utf8(&*aligned_seq1)?);
    println!("{}", from_utf8(&*diff_line)?);
    println!("{}", from_utf8(&*aligned_seq2)?);

    Ok(())
}
