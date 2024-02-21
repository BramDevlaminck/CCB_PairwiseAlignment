use std::cmp::{max, min};
use std::error::Error;
use std::str::from_utf8;
use clap::Parser;
use fasta_reader::read_fasta;

#[derive(Parser, Debug)]
#[clap(allow_negative_numbers = true)]
struct Args {
    /// The input file name
    #[clap(short, long)]
    filename: String,
    /// With of the band
    #[clap(short, long)]
    width: usize,
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
    let Args { filename, width, match_score, mismatch_score, gap_score } = args;
    let (seq1, seq2) = read_fasta(&filename)?;
    let seq1_chars = seq1.into_bytes();
    let seq2_chars = seq2.into_bytes();

    let mut matrix: Vec<Vec<Option<i32>>> = vec![vec![None; seq1_chars.len() + 1]; seq2_chars.len() + 1];

    // initialize the first column and first row
    matrix[0] = (0..=seq1_chars.len()).map(|i| Some(i as i32 * gap_score)).collect();
    for (index, init_value) in (0..=seq2_chars.len()).enumerate() {
        matrix[index][0] = Some(init_value as i32 * gap_score);
    }

    // fill in the matrix
    for row in 1..=seq2_chars.len() {
        let min_col = max(1, row as i32 - width as i32) as usize;
        let max_col = min(seq1_chars.len(), row + width);
        let current_seq2_char = seq2_chars[row - 1];
        for col in min_col..=max_col {
            let mut diag_score = i32::MIN;
            if let Some(old_diag) = matrix[row - 1][col - 1] {
                diag_score = old_diag + if seq1_chars[col - 1] == current_seq2_char { match_score } else { mismatch_score };
            }
            let mut col_score = i32::MIN;
            if let Some(old_col_score) = matrix[row][col - 1] {
                col_score = old_col_score + gap_score;
            }
            let mut row_score = i32::MIN;
            if let Some(old_row_score) = matrix[row - 1][col] {
                row_score = old_row_score + gap_score;
            }

            matrix[row][col] = Some(
                max(
                    max(diag_score, col_score),
                    row_score,
                )
            );
        }
    }

    let mut current_row = seq2_chars.len();
    let mut current_col = seq1_chars.len();
    println!("The score for optimal alignment is: {}", matrix[current_row][current_col].ok_or("Value in the right bottom corner was not calculated")?);

    // backtrack to find alignment
    let mut aligned_seq1: Vec<u8> = vec![];
    let mut aligned_seq2: Vec<u8> = vec![];
    let mut diff_line: Vec<u8> = vec![];

    let mut current_score = matrix[current_row][current_col].unwrap();
    while current_col != 0 && current_row != 0 {
        if matrix[current_row][current_col - 1].is_some() && current_col != 0 && current_score == matrix[current_row][current_col - 1].unwrap() + gap_score {
            aligned_seq1.insert(0, seq1_chars[current_col - 1]);
            diff_line.insert(0, b' ');
            aligned_seq2.insert(0, b'-');
            current_col -= 1;
        } else if matrix[current_row - 1][current_col].is_some() && current_row != 0 && current_score == matrix[current_row - 1][current_col].unwrap() + gap_score {
            aligned_seq1.insert(0, b'-');
            diff_line.insert(0, b' ');
            aligned_seq2.insert(0, seq2_chars[current_row - 1]);
            current_row -= 1;
        } else {
            // diagonal case
            aligned_seq1.insert(0, seq1_chars[current_col - 1]);
            aligned_seq2.insert(0, seq2_chars[current_row - 1]);
            let char_to_insert = if seq1_chars[current_col - 1] == seq2_chars[current_row - 1] { b'|' } else { b'*' };
            diff_line.insert(0, char_to_insert);
            current_row -= 1;
            current_col -= 1;
        }
        // update score
        current_score = matrix[current_row][current_col].unwrap();
    }

    println!();
    println!("Aligned sequences:");
    println!("{}", from_utf8(&*aligned_seq1)?);
    println!("{}", from_utf8(&*diff_line)?);
    println!("{}", from_utf8(&*aligned_seq2)?);

    Ok(())
}
