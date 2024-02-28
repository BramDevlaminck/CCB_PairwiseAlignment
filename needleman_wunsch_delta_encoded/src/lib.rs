pub fn construct_delta_matrices(seq1: &[u8], seq2: &[u8], match_score: i32, mismatch_score: i32, gap_score: i32) -> (Vec<Vec<Option<i32>>>, Vec<Vec<Option<i32>>>) {
    let mut delta_v: Vec<Vec<Option<i32>>> = vec![vec![None; seq1.len() + 1]; seq2.len() + 1];
    let mut delta_h: Vec<Vec<Option<i32>>> = vec![vec![None; seq1.len() + 1]; seq2.len() + 1];
    // initialize the first row and column
    delta_h[0] = (0..=seq1.len()).map(|_| Some(gap_score)).collect();
    for index in 0..=seq2.len() {
        delta_v[index][0] = Some(gap_score);
    }

    for row in 1..=seq2.len() {
        let current_seq2_char = seq2[row - 1];
        for col in 1..=seq1.len() {
            let cell_above = delta_h[row - 1][col].unwrap();
            let cell_left = delta_v[row][col - 1].unwrap();

            // fill in delta_v
            delta_v[row][col] = if seq1[col - 1] == current_seq2_char {
                Some(match_score - cell_above)
            } else if mismatch_score - gap_score >= cell_above && mismatch_score - gap_score >= cell_left { // mismatch
                Some(mismatch_score - cell_above)
            } else if cell_above >= mismatch_score - gap_score && cell_above >= cell_left { // indel from above
                Some(gap_score)
            } else { // indel from left
                assert!(cell_left >= mismatch_score - gap_score && cell_left >= cell_above);
                Some(cell_left + gap_score - cell_above)
            };

            // fill in delta_h
            delta_h[row][col] = if seq1[col - 1] == current_seq2_char {
                Some(match_score - cell_left)
            } else if mismatch_score - gap_score >= cell_above && mismatch_score - gap_score >= cell_left {
                Some(mismatch_score - cell_left)
            } else if cell_above >= mismatch_score - gap_score && cell_above >= cell_left { // indel from above
                Some(cell_above + gap_score - cell_left)
            } else { // indel from left
                assert!(cell_left >= mismatch_score - gap_score && cell_left >= cell_above);
                Some(gap_score)
            };

        }
    }

    (delta_h, delta_v)
}

#[cfg(test)]
mod tests {
    use std::error::Error;
    use fasta_reader::read_fasta;
    use needleman_wunsch::construct_matrix;
    use crate::construct_delta_matrices;

    #[test]
    fn test_delta_h_and_v_matrix_equal_to_s_matrix() -> Result<(), Box<dyn Error>> {
        let (seq1, seq2) = read_fasta("../tests/input.fasta")?;
        let seq1_chars = seq1.into_bytes();
        let seq2_chars = seq2.into_bytes();

        let match_score = 1;
        let mismatch_score = -1;
        let gap_score = -3;

        let (delta_h, delta_v) = construct_delta_matrices(&seq1_chars, &seq2_chars, match_score, mismatch_score, gap_score);

        let s_matrix = construct_matrix(&seq1_chars, &seq2_chars, match_score, mismatch_score, gap_score);

        // test delta_h matrix
        let mut s_matrix_from_h : Vec<Vec<i32>> = vec![vec![0; seq1_chars.len() + 1]; seq2_chars.len() + 1];
        for (index, init_value) in (0..=seq2_chars.len()).enumerate() {
            s_matrix_from_h[index][0] = init_value as i32 * gap_score;
        }
        for row in 0..=seq2_chars.len() {
            for col in 1..=seq1_chars.len() {
                s_matrix_from_h[row][col] = s_matrix_from_h[row][col-1] + delta_h[row][col].ok_or("Value should have been filled in")?;
            }
        }
        assert_eq!(s_matrix, s_matrix_from_h);

        let mut s_matrix_from_v : Vec<Vec<i32>> = vec![vec![0; seq1_chars.len() + 1]; seq2_chars.len() + 1];
        s_matrix_from_v[0] = (0..=seq1_chars.len()).map(|i| i as i32 * gap_score).collect();

        for row in 1..=seq2_chars.len() {
            for col in 0..=seq1_chars.len() {
                s_matrix_from_v[row][col] = s_matrix_from_v[row-1][col] + delta_v[row][col].ok_or("Value should have been filled in")?;
            }
        }
        assert_eq!(s_matrix, s_matrix_from_v);

        Ok(())
    }

}