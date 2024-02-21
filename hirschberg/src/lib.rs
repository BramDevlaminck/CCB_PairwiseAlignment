use std::cmp::max;
use std::error::Error;
use needleman_wunch::{backtrack_alignment, construct_matrix};


/// Calculate the needleman welsh score only using 2 rows
pub fn nw_score(seq1: &[u8], seq2: &[u8], match_score: i32, mismatch_score: i32, gap_score: i32, reversed: bool) -> Vec<i32> {
    let mut prev_row: Vec<i32> = (0..seq2.len() + 1).map(|i| i as i32 * gap_score).collect();
    for row in 1..=seq1.len() {
        let mut current_row = vec![row as i32 * gap_score];
        let current_seq1_char = if reversed {
            seq1[seq1.len() - row]
        } else {
            seq1[row - 1]
        };
        for col in 1..=seq2.len() {
            let current_seq2_char = if reversed {
                seq2[seq2.len() - col]
            } else {
                seq2[col - 1]
            };

            let diag_score = prev_row[col - 1] + if current_seq1_char == current_seq2_char { match_score } else { mismatch_score };
            current_row.push(
                max(
                    max(diag_score, current_row[col - 1] + gap_score),
                    prev_row[col] + gap_score,
                )
            );
        }
        prev_row = current_row;
    }
    prev_row
}

/// Execute the Hirschberg algorithm for global alignment on seq1 and seq2 using the provided scores
pub fn hirschberg(seq1: &[u8], seq2: &[u8], match_score: i32, mismatch_score: i32, gap_score: i32) -> Result<(Vec<u8>, Vec<u8>, Vec<u8>), Box<dyn Error>> {
    if seq1.len() <= 1 || seq2.len() <= 1 {
        let matrix = construct_matrix(&seq1, &seq2, match_score, mismatch_score, gap_score);
        return Ok(backtrack_alignment(&matrix, seq1, seq2, gap_score));
    }

    let xmid = seq1.len() / 2;

    let score_l = nw_score(&seq1[..xmid], &seq2, match_score, mismatch_score, gap_score, false);
    let mut score_r = nw_score(&seq1[xmid..], &seq2, match_score, mismatch_score, gap_score, true);
    score_r.reverse();

    let total_score = score_l.iter().zip(score_r).map(|(&l, r)| l + r).collect::<Vec<i32>>();
    let ymid = total_score.iter().enumerate().max_by_key(|(_, &key)| key).map(|(i, _)| i).ok_or("Score L and Score R are empty, no argmax can be found")?;

    let (mut alx1, mut diff1, mut aly1) = hirschberg(&seq1[..xmid], &seq2[..ymid], match_score, mismatch_score, gap_score)?;
    let (alx2, diff2, aly2) = hirschberg(&seq1[xmid..], &seq2[ymid..], match_score, mismatch_score, gap_score)?;

    alx1.extend(alx2);
    diff1.extend(diff2);
    aly1.extend(aly2);

    Ok((alx1, diff1, aly1))
}

#[cfg(test)]
mod tests {
    use std::error::Error;
    use fasta_reader::read_fasta;
    use crate::nw_score;

    #[test]
    fn test_nw_score() -> Result<(), Box<dyn Error>> {
        let (seq1, seq2) = read_fasta("../tests/input.fasta")?;
        let seq1_chars = seq1.into_bytes();
        let seq2_chars = seq2.into_bytes();

        let res = nw_score(&seq1_chars, &seq2_chars, 1, -1, -1, false);
        assert_eq!(res, vec![-7, -5, -3, -1, -2, -2, 0, 0]);
        Ok(())
    }


    #[test]
    fn test_nw_score_rev() -> Result<(), Box<dyn Error>> {
        let (seq1, seq2) = read_fasta("../tests/input.fasta")?;
        let seq1_chars = seq1.into_bytes();
        let seq2_chars = seq2.into_bytes();

        let res = nw_score(&seq1_chars, &seq2_chars, 1, -1, -1, true);
        assert_eq!(res, vec![-7, -7, -5, -3, -3, -1, -1, 0]);
        Ok(())
    }
}