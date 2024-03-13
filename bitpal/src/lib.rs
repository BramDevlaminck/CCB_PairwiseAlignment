use std::collections::{HashMap, HashSet};

use crate::bitpal_errors::InputTooLongError;

mod bitpal_errors;

/// The BitPAl algorithm implemented for scoring: M = 1, I = -1, G = -3
/// With the restriction that seq1 or seq2 needs to fit in 1 computer word (= seq1 or seq2 <= 64 characters)
pub fn bitpal(seq1: &Vec<u8>, seq2: &Vec<u8>) -> Result<i32, InputTooLongError> {
    // check validity of input
    let (ref horizontal_seq, vertical_seq) = if seq1.len() <= 64 {
        (seq1, seq2)
    } else if seq2.len() <= 64 {
        (seq2, seq1)
    } else {
        return Err(InputTooLongError);
    };

    // create a set containing the used alphabet
    let mut alphabet: HashSet<u8> = HashSet::new();
    alphabet.extend(horizontal_seq.iter());
    alphabet.extend(vertical_seq.iter());
    // build the needed match vectors
    let match_vectors = calculate_match_vectors(&horizontal_seq, &alphabet);

    // vector containing horizontal_seq.len() 1's
    let all_ones: u64 = 2_u64.wrapping_pow(horizontal_seq.len() as u32).wrapping_sub(1);
    let mut delta_h_pos4: u64 = 0;
    let mut delta_h_pos3: u64 = 0;
    let mut delta_h_pos2: u64 = 0;
    let mut delta_h_pos1: u64 = 0;
    let mut delta_h_0: u64 = 0;
    let mut delta_h_neg1: u64 = 0;
    let mut delta_h_neg2: u64 = 0;
    let mut delta_h_neg3: u64 = all_ones; // initialize first delta H with all 1's

    for &character in vertical_seq.iter() {
        let current_match_vector = match_vectors.get(&character).unwrap();
        let not_match = !current_match_vector;
        // calculate max value
        let init_pos4 = current_match_vector & delta_h_neg3;
        let delta_v_pos4_shift = ((init_pos4.wrapping_add(delta_h_neg3)) ^ delta_h_neg3) ^ init_pos4;

        // calculate rest of delta v high from high to low
        let remain_delta_h_neg3 = delta_h_neg3 ^ (delta_v_pos4_shift >> 1);
        let delta_v_pos4_shift_or_match = delta_v_pos4_shift | current_match_vector;
        let init_pos3s = delta_h_neg2 & delta_v_pos4_shift_or_match;
        let delta_v_pos3_shift = ((init_pos3s << 1).wrapping_add(remain_delta_h_neg3)) ^ remain_delta_h_neg3;
        let delta_v_pos3_shift_not_match = delta_v_pos3_shift & not_match;

        // calculate min + 1 to mid
        let delta_v_not_4_to_3_shift_or_match = !(delta_v_pos4_shift_or_match | delta_v_pos3_shift);
        let delta_v_pos2_shift = ((delta_v_pos4_shift_or_match & delta_h_neg1) | (delta_v_pos3_shift_not_match & delta_h_neg2) | (delta_v_not_4_to_3_shift_or_match & delta_h_neg3)) << 1;
        let delta_v_pos1_shift = ((delta_v_pos4_shift_or_match & delta_h_0) | (delta_v_pos3_shift_not_match & delta_h_neg1) | (delta_v_not_4_to_3_shift_or_match & delta_h_neg2)) << 1;
        let delta_v_0_shift = ((delta_v_pos4_shift_or_match & delta_h_pos1) | (delta_v_pos3_shift_not_match & delta_h_0) | (delta_v_not_4_to_3_shift_or_match & delta_h_neg1)) << 1;
        let delta_v_neg1_shift = ((delta_v_pos4_shift_or_match & delta_h_pos2) | (delta_v_pos3_shift_not_match & delta_h_pos1) | (delta_v_not_4_to_3_shift_or_match & delta_h_0)) << 1;
        let delta_v_neg2_shift = ((delta_v_pos4_shift_or_match & delta_h_pos3) | (delta_v_pos3_shift_not_match & delta_h_pos2) | (delta_v_not_4_to_3_shift_or_match & delta_h_pos1)) << 1;

        // fill in remaining min values
        let delta_v_neg3_shift = all_ones ^ (delta_v_pos4_shift | delta_v_pos3_shift | delta_v_pos2_shift | delta_v_pos1_shift | delta_v_0_shift | delta_v_neg1_shift | delta_v_neg2_shift);

        // prepare the delta H vectors for calculation
        delta_h_pos4 |= current_match_vector; // add matches to the max vector
        delta_h_pos2 = (delta_h_pos2 | delta_h_pos1 | delta_h_0 | delta_h_neg1 | delta_h_neg2 | delta_h_neg3) & not_match; // low to mid and remove match
        delta_h_pos3 = delta_h_pos3 & not_match; // remove match

        // calculate the new delta h vectors (from low to high since we use delta_h_pos2, delta_h_pos3 and delta_h_pos4 in the lowest values, so we can only modify them in the end)
        delta_h_neg2 = (delta_h_pos4 & delta_v_pos3_shift) | (delta_h_pos3 & delta_v_pos2_shift) | (delta_h_pos2 & delta_v_pos1_shift);
        delta_h_neg1 = (delta_h_pos4 & delta_v_pos2_shift) | (delta_h_pos3 & delta_v_pos1_shift) | (delta_h_pos2 & delta_v_0_shift);
        delta_h_0 = (delta_h_pos4 & delta_v_pos1_shift) | (delta_h_pos3 & delta_v_0_shift) | (delta_h_pos2 & delta_v_neg1_shift);
        delta_h_pos1 = (delta_h_pos4 & delta_v_0_shift) | (delta_h_pos3 & delta_v_neg1_shift) | (delta_h_pos2 & delta_v_neg2_shift);
        delta_h_pos2 = (delta_h_pos4 & delta_v_neg1_shift) | (delta_h_pos3 & delta_v_neg2_shift) | (delta_h_pos2 & delta_v_neg3_shift);
        delta_h_pos3 = (delta_h_pos4 & delta_v_neg2_shift) | (delta_h_pos3 & delta_v_neg3_shift);
        delta_h_pos4 &= delta_v_neg3_shift;

        // fill in remaining min values
        delta_h_neg3 = all_ones ^ (delta_h_pos4 | delta_h_pos3 | delta_h_pos2 | delta_h_pos1 | delta_h_0 | delta_h_neg1 | delta_h_neg2);
    }

    let score_mapping = HashMap::from([
        (-3, delta_h_neg3),
        (-2, delta_h_neg2),
        (-1, delta_h_neg1),
        (1, delta_h_pos1),
        (2, delta_h_pos2),
        (3, delta_h_pos3),
        (4, delta_h_pos4)
    ]);

    Ok(gather_score(score_mapping, vertical_seq.len(), -3))
}

#[allow(unused)]
fn format_binary(vector: u64) -> String {
    format!("{:064b}", vector).chars().rev().collect::<String>()
}

#[allow(unused)]
fn print_binary(vector: u64) {
    println!("{}", format_binary(vector));
}

/// Calculate a hashmap with a bitvector representing the locations where a character occurs in seq.
fn calculate_match_vectors(seq: &Vec<u8>, alphabet: &HashSet<u8>) -> HashMap<u8, u64> {
    let mut match_vectors: HashMap<u8, u64> = HashMap::new();
    for &letter in alphabet {
        match_vectors.insert(letter, 0);
    }

    for (i, &character) in seq.iter().enumerate() {
        match_vectors.entry(character).and_modify(|vector| *vector |= 1 << i);
    }

    match_vectors
}

/// Calculate the score using the bitvectors provided in `delta_h_map`.
/// The key is the value i, that the bits represent, while the value is the bitvector
fn gather_score(delta_h_map: HashMap<i32, u64>, vertical_sequence_length: usize, gap_score: i32) -> i32 {
    let mut score = gap_score * vertical_sequence_length as i32;
    for (value, bitrow) in delta_h_map {
        score += bitrow.count_ones() as i32 * value;
    }
    score
}

#[cfg(test)]
mod tests {
    use std::collections::{HashMap, HashSet};
    use std::error::Error;
    use std::ops::Range;

    use rand::Rng;
    use rand::rngs::ThreadRng;

    use fasta_reader::read_fasta;
    use needleman_wunsch::construct_matrix;

    use crate::{bitpal, calculate_match_vectors, gather_score};
    use crate::bitpal_errors::InputTooLongError;

    #[test]
    fn test_match_vector() {
        let seq = "ACACGTA".to_string();
        let alphabet = HashSet::from([b'A', b'C', b'G', b'T']);
        let res = calculate_match_vectors(&seq.into_bytes(), &alphabet);
        let expected = HashMap::from([
            (b'A', 69),
            (b'C', 10),
            (b'G', 16),
            (b'T', 32),
        ]);
        assert_eq!(res, expected);
    }

    #[test]
    fn test_gather_score() {
        let last_row_res = HashMap::from([
            (-3, 4),
            (-2, 6),
            (-1, 0),
            (1, 0),
            (2, 10),
            (3, 69),
            (4, 3)
        ]);

        let res = gather_score(last_row_res, 5, -3);

        let expected = 5 * -3 + (1 * -3 + 2 * -2 + 2 * 2 + 3 * 3 + 2 * 4);
        assert_eq!(res, expected);
    }

    #[test]
    fn test_bitpal() -> Result<(), Box<dyn Error>> {
        let (seq1, seq2) = read_fasta("../tests/dna_input.fasta")?;
        let seq1_chars = seq1.into_bytes();
        let seq2_chars = seq2.into_bytes();

        let matrix = construct_matrix(&seq1_chars, &seq2_chars, 1, -1, -3);
        let nw_score = matrix[seq2_chars.len()][seq1_chars.len()];
        let bitpal = bitpal(&seq1_chars, &seq2_chars)?;

        assert_eq!(bitpal, nw_score);

        Ok(())
    }

    #[test]
    fn test_bitpal_minimal() -> Result<(), Box<dyn Error>> {
        let seq1 = "A".to_string().into_bytes();
        let seq2 = "TG".to_string().into_bytes();
        let matrix = construct_matrix(&seq1, &seq2, 1, -1, -3);
        let nw_score = matrix[seq2.len()][seq1.len()];
        let bitpal = bitpal(&seq1, &seq2)?;

        assert_eq!(bitpal, nw_score);

        Ok(())
    }

    /// Helper function to fuzz test BitPAl implementation where a random sequence with a length in the size_range is generated using the provided alphabet
    fn generate_sequence_with_alphabet(alphabet: &Vec<u8>, rng: &mut ThreadRng, size_range: Range<usize>) -> Vec<u8> {
        let seq_len = rng.gen_range(size_range);
        let mut seq = vec![];
        for _ in 0..seq_len {
            seq.push(alphabet[rng.gen_range(0..alphabet.len())]);
        }

        seq
    }

    #[test]
    fn test_bitpal_fuzzing() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::thread_rng();
        let valid_letters = vec![b'A', b'C', b'G', b'T'];
        for _ in 0..5000 {
            let seq1 = generate_sequence_with_alphabet(&valid_letters, &mut rng, 1..65); // this one is set horizontal, so make sure it is not longer than 1 word
            let seq2 = generate_sequence_with_alphabet(&valid_letters, &mut rng, 1..513);
            let matrix = construct_matrix(&seq1, &seq2, 1, -1, -3);
            let nw_score = matrix[seq2.len()][seq1.len()];
            let bitpal = bitpal(&seq1, &seq2)?;

            assert_eq!(bitpal, nw_score);
        }

        Ok(())
    }

    #[test]
    fn test_bitpal_size_check() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::thread_rng();
        let valid_letters = vec![b'A', b'C', b'G', b'T'];
        // make both sequences larger than 65
        let seq1 = generate_sequence_with_alphabet(&valid_letters, &mut rng, 65..66);
        let seq2 = generate_sequence_with_alphabet(&valid_letters, &mut rng, 65..66);
        assert!(bitpal(&seq1, &seq2).is_err_and(|e| e == InputTooLongError));

        Ok(())
    }

    #[test]
    fn test_bitpal_seq2_short() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::thread_rng();
        let valid_letters = vec![b'A', b'C', b'G', b'T'];
        // make seq2 the sequence that should be horizontal
        let seq1 = generate_sequence_with_alphabet(&valid_letters, &mut rng, 65..66);
        let seq2 = generate_sequence_with_alphabet(&valid_letters, &mut rng, 1..65);
        let matrix = construct_matrix(&seq1, &seq2, 1, -1, -3);
        let nw_score = matrix[seq2.len()][seq1.len()];
        let bitpal = bitpal(&seq1, &seq2)?;
        assert_eq!(bitpal, nw_score);

        Ok(())
    }

    #[test]
    fn test_bitpal_seq1_empty() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::thread_rng();
        let valid_letters = vec![b'A', b'C', b'G', b'T'];
        // make seq2 the sequence that should be horizontal
        let seq1 = vec![];
        let seq2 = generate_sequence_with_alphabet(&valid_letters, &mut rng, 1..65);
        let bitpal = bitpal(&seq1, &seq2)?;
        assert_eq!(bitpal, -3 * seq2.len() as i32);

        Ok(())
    }

    #[test]
    fn test_bitpal_seq2_empty() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::thread_rng();
        let valid_letters = vec![b'A', b'C', b'G', b'T'];
        // make seq2 the sequence that should be horizontal
        let seq1 = generate_sequence_with_alphabet(&valid_letters, &mut rng, 1..65);
        let seq2 = vec![];
        let bitpal = bitpal(&seq1, &seq2)?;
        assert_eq!(bitpal, -3 * seq1.len() as i32);

        Ok(())
    }

    #[test]
    fn test_bitpal_seq1_empty_seq2_empty() -> Result<(), Box<dyn Error>> {
        // make seq2 the sequence that should be horizontal
        let seq1 = vec![];
        let seq2 = vec![];
        let bitpal = bitpal(&seq1, &seq2)?;
        assert_eq!(bitpal, 0);

        Ok(())
    }

    #[test]
    fn test_bitpal_seq_len64() -> Result<(), Box<dyn Error>> {
        let mut rng = rand::thread_rng();
        // make seq2 the sequence that should be horizontal
        let valid_letters = vec![b'A', b'C', b'G', b'T'];
        let seq1 = generate_sequence_with_alphabet(&valid_letters, &mut rng, 64..65);
        let seq2 = generate_sequence_with_alphabet(&valid_letters, &mut rng, 64..65);
        let matrix = construct_matrix(&seq1, &seq2, 1, -1, -3);
        let nw_score = matrix[seq2.len()][seq1.len()];
        let bitpal = bitpal(&seq1, &seq2)?;
        assert_eq!(bitpal, nw_score);

        Ok(())
    }
}