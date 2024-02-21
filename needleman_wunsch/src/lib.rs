use std::cmp::max;

pub fn construct_matrix(seq1: &[u8], seq2: &[u8], match_score: i32, mismatch_score: i32, gap_score: i32) -> Vec<Vec<i32>> {
    let mut matrix: Vec<Vec<i32>> = vec![vec![0; seq1.len() + 1]; seq2.len() + 1];

    // initialize the first column and first row
    matrix[0] = (0..=seq1.len()).map(|i| i as i32 * gap_score).collect();
    for (index, init_value) in (0..seq2.len() + 1).enumerate() {
        matrix[index][0] = init_value as i32 * gap_score;
    }

    // fill in the matrix
    for row in 1..=seq2.len() {
        let current_seq2_char = seq2[row - 1];
        for col in 1..=seq1.len() {
            let diag_score = matrix[row - 1][col - 1] + if seq1[col - 1] == current_seq2_char { match_score } else { mismatch_score };
            matrix[row][col] = max(
                max(diag_score, matrix[row][col - 1] + gap_score),
                matrix[row - 1][col] + gap_score,
            );
        }
    }

    matrix
}

pub fn backtrack_alignment(matrix: &Vec<Vec<i32>>, seq1: &[u8], seq2: &[u8], gap_score: i32) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let mut current_row = seq2.len();
    let mut current_col = seq1.len();

    let mut aligned_seq1: Vec<u8> = vec![];
    let mut aligned_seq2: Vec<u8> = vec![];
    let mut diff_line: Vec<u8> = vec![];

    while current_col != 0 || current_row != 0 {
        if current_col != 0 && matrix[current_row][current_col] == matrix[current_row][current_col - 1] + gap_score {
            aligned_seq1.insert(0, seq1[current_col - 1]);
            diff_line.insert(0, b' ');
            aligned_seq2.insert(0, b'-');
            current_col -= 1;
            continue;
        }

        if current_row != 0 && matrix[current_row][current_col] == matrix[current_row - 1][current_col] + gap_score {
            aligned_seq1.insert(0, b'-');
            diff_line.insert(0, b' ');
            aligned_seq2.insert(0, seq2[current_row - 1]);
            current_row -= 1;
            continue;
        }

        // diagonal case
        aligned_seq1.insert(0, seq1[current_col - 1]);
        aligned_seq2.insert(0, seq2[current_row - 1]);
        let char_to_insert = if seq1[current_col - 1] == seq2[current_row - 1] { b'|' } else { b'*' };
        diff_line.insert(0, char_to_insert);
        current_row -= 1;
        current_col -= 1;
    }

    (aligned_seq1, diff_line, aligned_seq2)
}


#[cfg(test)]
mod tests {
    use std::error::Error;
    use fasta_reader::read_fasta;
    use crate::{backtrack_alignment, construct_matrix};

    #[test]
    fn test_matrix() -> Result<(), Box<dyn Error>> {
        let (seq1, seq2) = read_fasta("../tests/input.fasta")?;
        let seq1_chars = seq1.into_bytes();
        let seq2_chars = seq2.into_bytes();

        let res = construct_matrix(&seq1_chars, &seq2_chars, 1, -1, -1);
        assert_eq!(res[res.len() - 1], vec![-7, -5, -3, -1, -1, -1, 0, 0]);
        Ok(())
    }


    #[test]
    fn test_backtrack() -> Result<(), Box<dyn Error>> {
        let (seq1, seq2) = read_fasta("../tests/input.fasta")?;
        let seq1_chars = seq1.into_bytes();
        let seq2_chars = seq2.into_bytes();

        let matrix = construct_matrix(&seq1_chars, &seq2_chars, 1, -1, -1);
        let res = backtrack_alignment(&matrix, &seq1_chars, &seq2_chars, -1);
        assert_eq!(res, ("G-ATTACA".to_string().into_bytes(), "| ||* |*".to_string().into_bytes(), "GCATG-CU".to_string().into_bytes()));
        Ok(())
    }

    #[test]
    fn test_backtrack_seq2_empty() -> Result<(), Box<dyn Error>> {
        let seq1 = "test".to_string();
        let seq2 = String::new();
        let seq1_chars = seq1.into_bytes();
        let seq2_chars = seq2.into_bytes();

        let matrix = construct_matrix(&seq1_chars, &seq2_chars, 1, -1, -1);
        let res = backtrack_alignment(&matrix, &seq1_chars, &seq2_chars, -1);
        assert_eq!(res, ("test".to_string().into_bytes(), "    ".to_string().into_bytes(), "----".to_string().into_bytes()));
        Ok(())
    }

    #[test]
    fn test_backtrack_seq1_empty() -> Result<(), Box<dyn Error>> {
        let seq1 = String::new();
        let seq2 = "test".to_string();
        let seq1_chars = seq1.into_bytes();
        let seq2_chars = seq2.into_bytes();

        let matrix = construct_matrix(&seq1_chars, &seq2_chars, 1, -1, -1);
        let res = backtrack_alignment(&matrix, &seq1_chars, &seq2_chars, -1);
        assert_eq!(res, ("----".to_string().into_bytes(), "    ".to_string().into_bytes(), "test".to_string().into_bytes()));
        Ok(())
    }
}