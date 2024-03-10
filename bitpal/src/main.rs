use std::collections::HashMap;
use std::error::Error;
use std::hash::Hash;
use clap::Parser;
use fasta_reader::read_fasta;

#[derive(Parser, Debug)]
#[clap(allow_negative_numbers = true)]
struct Args {
    /// The input file name
    #[clap(short, long)]
    filename: String,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let Args { filename } = args;
    let (seq1, seq2) = read_fasta(&filename)?;
    let seq1_chars = seq1.into_bytes();
    let seq2_chars = seq2.into_bytes();

    // TODO: some check that at least 1 should be max 64 characters

    let result_score = bitpal(seq1_chars, seq2_chars);
    println!("{}", result_score);

    Ok(())
}

fn bitpal(seq1: Vec<u8>, seq2: Vec<u8>) -> i32 {
    let all_ones: u64 = 2_u64.pow(seq1.len() as u32) - 1;
    let match_vectors = calculate_match_vectors(&seq1);

    let mut delta_h_pos4: u64 = 0;
    let mut delta_h_pos3: u64 = 0;
    let mut delta_h_pos2: u64 = 0;
    let mut delta_h_pos1: u64 = 0;
    let mut delta_h_0: u64 = 0;
    let mut delta_h_neg1: u64 = 0;
    let mut delta_h_neg2: u64 = 0;
    let mut delta_h_neg3: u64 = all_ones; // initialize first delta H with all 1's

    for &character in seq2.iter() {
        let current_match_vector = match_vectors.get(&(character as char)).unwrap();
        let current_match_vector_format = format_binary(*current_match_vector);
        let not_match = !current_match_vector;
        let not_match_format = format_binary(not_match);
        // calculate max value
        let init_pos4 = current_match_vector & delta_h_neg3;
        let delta_v_pos4_shift = ((init_pos4 + delta_h_neg3) ^ delta_h_neg3) ^ init_pos4;
        let delta_v_pos4 = delta_v_pos4_shift >> 1;

        // calculate rest of delta v high from high to low
        let remain_delta_h_neg3 = delta_h_neg3 ^ delta_v_pos4;
        let delta_v_pos4_shift_or_match = delta_v_pos4_shift | current_match_vector;
        let init_pos3s = delta_h_neg2 & delta_v_pos4_shift_or_match;
        let delta_v_pos3_shift = ((init_pos3s << 1) + remain_delta_h_neg3) ^ remain_delta_h_neg3;
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

        // fill in the delta H values
        delta_h_pos4 |= current_match_vector; // add matches to the max vector
        let delta_h_pos2_no_match = (delta_h_pos2 | delta_h_pos1 | delta_h_0 | delta_h_neg1 | delta_h_neg2 | delta_h_neg3) & not_match; // low to mid and remove match
        let delta_h_pos3_no_match = delta_h_pos3 & not_match; // remove match

        let delta_h_pos4_new = delta_h_pos4 & delta_v_neg3_shift;
        delta_h_pos3 = (delta_h_pos4 & delta_v_neg2_shift) | (delta_h_pos3_no_match & delta_v_neg3_shift);
        delta_h_pos2 = (delta_h_pos4 & delta_v_neg1_shift) | (delta_h_pos3_no_match & delta_v_neg2_shift) | (delta_h_pos2_no_match & delta_v_neg3_shift);
        delta_h_pos1 = (delta_h_pos4 & delta_v_0_shift) | (delta_h_pos3_no_match & delta_v_neg1_shift) | (delta_h_pos2_no_match & delta_v_neg2_shift);
        delta_h_0 = (delta_h_pos4 & delta_v_pos1_shift) | (delta_h_pos3_no_match & delta_v_0_shift) | (delta_h_pos2_no_match & delta_v_neg1_shift);
        delta_h_neg1 = (delta_h_pos4 & delta_v_pos2_shift) | (delta_h_pos3_no_match & delta_v_pos1_shift) | (delta_h_pos2_no_match & delta_v_0_shift);
        delta_h_neg2 = (delta_h_pos4 & delta_v_pos3_shift) | (delta_h_pos3_no_match & delta_v_pos2_shift) | (delta_h_pos2_no_match & delta_v_pos1_shift);

        delta_h_pos4 = delta_h_pos4_new;

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

    gather_score(score_mapping, seq2.len(), -3)
}

fn format_binary(vector: u64) -> String {
    format!("{:064b}", vector).chars().rev().collect::<String>()
}

fn print_binary(vector: u64) {
    println!("{}", format_binary(vector));
}

/// Calculate a hashmap with a bitvector representing the locations where a character occurs in seq.
fn calculate_match_vectors(seq: &Vec<u8>) -> HashMap<char, u64> {
    // put seq 1 horizontally, so calculate the match vectors for it
    let mut match_vectors: HashMap<char, u64> = HashMap::from([
        ('A', 0),
        ('C', 0),
        ('G', 0),
        ('T', 0),
    ]);

    for (i, &character) in seq.iter().enumerate() {
        match_vectors.entry(character as char).and_modify(|vector| *vector |= 1 << i);
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
    use std::collections::HashMap;
    use std::error::Error;
    use fasta_reader::read_fasta;
    use needleman_wunsch::construct_matrix;
    use crate::{bitpal, calculate_match_vectors, gather_score};
    use rand::Rng;
    use rand::rngs::ThreadRng;


    #[test]
    fn test_match_vector() {
        let seq = "ACACGTA".to_string();
        let res = calculate_match_vectors(&seq.into_bytes());
        let expected = HashMap::from([
            ('A', 69),
            ('C', 10),
            ('G', 16),
            ('T', 32),
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
        let bitpal = bitpal(seq1_chars, seq2_chars);

        assert_eq!(bitpal, nw_score);

        Ok(())
    }

    fn generate_sequence_with_alphabet(alphabet: &Vec<u8>, rng: &mut ThreadRng) -> Vec<u8> {
        let seq_len = rng.gen_range(1..64);
        let mut seq = vec![];
        for _ in 0..seq_len {
            seq.push(alphabet[rng.gen_range(0..alphabet.len())]);
        }

        seq
    }

    #[test]
    fn test_bitpal_fuzzing() {
        let mut rng = rand::thread_rng();
        let valid_letters = vec![b'A', b'C', b'G', b'T'];
        for _ in 0..1000 {
            let seq1 = generate_sequence_with_alphabet(&valid_letters, &mut rng);
            let seq2 = generate_sequence_with_alphabet(&valid_letters, &mut rng);
            let matrix = construct_matrix(&seq1, &seq2, 1, -1, -3);
            let nw_score = matrix[seq2.len()][seq1.len()];
            let bitpal = bitpal(seq1, seq2);

            assert_eq!(bitpal, nw_score);
        }
    }

    #[test]
    fn test_bitpal_minimal() {
        let seq1 = "A".to_string().into_bytes();
        let seq2 = "TG".to_string().into_bytes();
        let matrix = construct_matrix(&seq1, &seq2, 1, -1, -3);
        let nw_score = matrix[seq2.len()][seq1.len()];
        let bitpal = bitpal(seq1, seq2);

        assert_eq!(bitpal, nw_score);
    }
}