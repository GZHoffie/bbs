use crate::types::*;
use std::collections::HashSet;
use std::collections::HashMap;


#[derive(Debug, Clone, PartialEq, Eq, Hash, Copy)]
pub struct KMer {
    k: usize,

    // represent k-mer using 3 u32 objects.
    seq1: u32,
    seq2: u32,
    seq3: u32,
}

impl KMer {
    pub fn new(k: usize) -> Self {
        if k > 64 {
            panic!("K-mer length too long to fit in 3 u32 objects");
        }


        KMer {
            k,
            seq1: 0,
            seq2: 0,
            seq3: 0,
        }
    }

    pub fn from_str(k: usize, kmer_str: &str) -> Self {
        if kmer_str.len() > k {
            panic!("K-mer string length does not match k");
        }

        let mut kmer = KMer::new(k);

        for (i, b) in kmer_str.bytes().enumerate() {
            let nuc = BYTE_TO_SEQ[b as usize] as u32;
            
            let nuc1 = nuc & 1;
            let nuc2 = (nuc >> 1) & 1;
            let nuc3 = (nuc >> 2) & 1;

            kmer.seq1 = (kmer.seq1 << 1) | nuc1;
            kmer.seq2 = (kmer.seq2 << 1) | nuc2;
            kmer.seq3 = (kmer.seq3 << 1) | nuc3;
        }

        kmer
    }

    /**
     * Return the last nucleotide of the k-mer as a u32, where A=0, C=1, G=2, T=3.
     */
    pub fn last_nuc_as_u32(&self) -> u32 {
        let nuc1 = self.seq1 & 1;
        let nuc2 = self.seq2 & 1;
        let nuc3 = self.seq3 & 1;
        
        (nuc3 << 2) | (nuc2 << 1) | nuc1
    }


    pub fn first_nuc_as_u32(&self) -> u32 {
        let shift = self.k - 1;
        let nuc1 = (self.seq1 >> shift) & 1;
        let nuc2 = (self.seq2 >> shift) & 1;
        let nuc3 = (self.seq3 >> shift) & 1;
        
        (nuc3 << 2) | (nuc2 << 1) | nuc1
    }


    /**
     * KMer objects that can be the next k-mer in the de Bruijn graph, 
     * by shifting the current k-mer and adding a new nucleotide at the end.
     */
    pub fn next_kmers(&self) -> Vec<KMer> {
        let mut next_kmers = Vec::new();

        // A, C, G, T
        for nuc in 0..4 {
            let nuc1 = nuc & 1;
            let nuc2 = (nuc >> 1) & 1;
            let nuc3 = (nuc >> 2) & 1;

            let next_seq1 = ((self.seq1 << 1) | nuc1) & ((1 << self.k) - 1);
            let next_seq2 = ((self.seq2 << 1) | nuc2) & ((1 << self.k) - 1);
            let next_seq3 = ((self.seq3 << 1) | nuc3) & ((1 << self.k) - 1);

            next_kmers.push(KMer {
                k: self.k,
                seq1: next_seq1,
                seq2: next_seq2,
                seq3: next_seq3,
            });
        }

        next_kmers
    }

    /**
     * Util function for the KMC class.
     * Return the next k-mers of self, but
     * -> if self == ending_kmer, then return an empty vector, which represents the end of the sequence.
     * -> if self ends with '$', then return a vector containing only the k-mer with '$' added at the end, which represents the end of the sequence.
     * -> Otherwise, return next_kmers() plus the k-mer with '$' added at the end.
     */
    pub fn next_kmers_at_end(&self, ending_kmer: KMer) -> Vec<KMer> {
        if *self == ending_kmer {
            return Vec::new();
        } else if self.last_nuc_as_u32() == END_MARKER as u32 {
            return vec![self.next_kmer_end()];
        }

        let mut next_kmers = self.next_kmers();

        // add the k-mer with '$' added at the end
        let next_kmer_end = self.next_kmer_end();
        next_kmers.push(next_kmer_end);

        next_kmers
    }



    /**
     * Neighboring k-mers that is same as self except for the last base, substituted to
     * A, C, G, T.
     */
    pub fn neighbor_kmers_last_base(&self) -> Vec<KMer> {
        let mut next_kmers = Vec::new();

        // A, C, G, T
        for nuc in 0..4 {
            let nuc1 = nuc & 1;
            let nuc2 = (nuc >> 1) & 1;
            let nuc3 = (nuc >> 2) & 1;

            let next_seq1 = ((self.seq1 << 1) | nuc1) & ((1 << self.k) - 1);
            let next_seq2 = ((self.seq2 << 1) | nuc2) & ((1 << self.k) - 1);
            let next_seq3 = ((self.seq3 << 1) | nuc3) & ((1 << self.k) - 1);

            next_kmers.push(KMer {
                k: self.k,
                seq1: next_seq1,
                seq2: next_seq2,
                seq3: next_seq3,
            });
        }

        next_kmers
    }

    /**
     * Neighboring k-mers that is same as self except for the first base, substituted to
     * A, C, G, T.
     */
    pub fn neighbor_kmers_first_base(&self) -> Vec<KMer> {
        let mut prev_kmers = Vec::new();

        // A, C, G, T
        for nuc in 0..4 {
            let nuc1 = nuc & 1;
            let nuc2 = (nuc >> 1) & 1;
            let nuc3 = (nuc >> 2) & 1;

            let prev_seq1 = ((self.seq1 >> 1) | (nuc1 << (self.k - 1))) & ((1 << self.k) - 1);
            let prev_seq2 = ((self.seq2 >> 1) | (nuc2 << (self.k - 1))) & ((1 << self.k) - 1);
            let prev_seq3 = ((self.seq3 >> 1) | (nuc3 << (self.k - 1))) & ((1 << self.k) - 1);

            prev_kmers.push(KMer {
                k: self.k,
                seq1: prev_seq1,
                seq2: prev_seq2,
                seq3: prev_seq3,
            });
        }

        prev_kmers
    }

    /**
     * A KMer object that can be the next k-mer in the de Bruijn graph, 
     * by shifting the current k-mer and adding an end marker ('$') at the end.
     */
    pub fn next_kmer_end(&self) -> KMer {
        let nuc1 = (END_MARKER & 1) as u32;
        let nuc2 = ((END_MARKER >> 1) & 1) as u32;
        let nuc3 = ((END_MARKER >> 2) & 1) as u32;

        let next_seq1 = ((self.seq1 << 1) | nuc1) & ((1 << self.k) - 1);
        let next_seq2 = ((self.seq2 << 1) | nuc2) & ((1 << self.k) - 1);
        let next_seq3 = ((self.seq3 << 1) | nuc3) & ((1 << self.k) - 1);

        KMer {
            k: self.k,
            seq1: next_seq1,
            seq2: next_seq2,
            seq3: next_seq3,
        }
    }


    fn _update_with_nuc(&mut self, nuc: u8) {
        let nuc_f = BYTE_TO_SEQ[nuc as usize] as u32;

        let nuc1 = nuc_f & 1;
        let nuc2 = (nuc_f >> 1) & 1;
        let nuc3 = (nuc_f >> 2) & 1;

        self.seq1 = ((self.seq1 << 1) | nuc1) & ((1 << self.k) - 1);
        self.seq2 = ((self.seq2 << 1) | nuc2) & ((1 << self.k) - 1);
        self.seq3 = ((self.seq3 << 1) | nuc3) & ((1 << self.k) - 1);
    }

    pub fn prev_kmers(&self) -> Vec<KMer> {
        let mut prev_kmers = Vec::new();

        // A, C, G, T
        for nuc in 0..4 {
            let nuc1 = nuc & 1;
            let nuc2 = (nuc >> 1) & 1;
            let nuc3 = (nuc >> 2) & 1;

            let prev_seq1 = ((self.seq1 >> 1) | (nuc1 << (self.k - 1))) & ((1 << self.k) - 1);
            let prev_seq2 = ((self.seq2 >> 1) | (nuc2 << (self.k - 1))) & ((1 << self.k) - 1);
            let prev_seq3 = ((self.seq3 >> 1) | (nuc3 << (self.k - 1))) & ((1 << self.k) - 1);

            prev_kmers.push(KMer {
                k: self.k,
                seq1: prev_seq1,
                seq2: prev_seq2,
                seq3: prev_seq3,
            });
        }

        prev_kmers
    }

    /**
     * Same logic as next_kmers_at_end
     */
    pub fn prev_kmers_at_end(&self, ending_kmer: KMer) -> Vec<KMer> {
        if *self == ending_kmer {
            return Vec::new();
        } else if self.first_nuc_as_u32() == END_MARKER as u32 {
            return vec![self.prev_kmer_end()];
        }

        let mut prev_kmers = self.prev_kmers();

        prev_kmers.push(self.prev_kmer_end());

        prev_kmers
    }

    pub fn prev_kmer_end(&self) -> KMer {
        let nuc1 = (END_MARKER & 1) as u32;
        let nuc2 = ((END_MARKER >> 1) & 1) as u32;
        let nuc3 = ((END_MARKER >> 2) & 1) as u32;

        let prev_seq1 = ((self.seq1 >> 1) | (nuc1 << (self.k - 1))) & ((1 << self.k) - 1);
        let prev_seq2 = ((self.seq2 >> 1) | (nuc2 << (self.k - 1))) & ((1 << self.k) - 1);
        let prev_seq3 = ((self.seq3 >> 1) | (nuc3 << (self.k - 1))) & ((1 << self.k) - 1);

        KMer {
            k: self.k,
            seq1: prev_seq1,
            seq2: prev_seq2,
            seq3: prev_seq3,
        }
    }

    /**
     * Convert the KMer object back to a string representation.
     */
    pub fn to_string(&self) -> String {
        let mut kmer_str = String::new();
        for i in (0..self.k).rev() {
            let nuc1 = (self.seq1 >> i) & 1;
            let nuc2 = (self.seq2 >> i) & 1;
            let nuc3 = (self.seq3 >> i) & 1;
            let nuc = (nuc3 << 2) | (nuc2 << 1) | nuc1;
            kmer_str.push(SEQ_TO_BYTE[nuc as usize] as char);

        }
        kmer_str
    }
}



/**
 * @brief Compute the k-mer vector for a given string. Only consider
 * the forward strand.
 * 
 * @param seq The string to compute the k-mer vector for
 * @param k The length of the k-mers
 * @returns A vector of k-mers
 */
pub fn seq_to_kmer_vec(
    seq: &String,
    k: usize,
    append_end: bool,
) -> Vec<KMer> {
    // If the string is shorter than the k-mer length, return
    if seq.len() < k {
        return Vec::new();
    }

    // Append k '$' characters at the front and back of the string to handle edge cases
    let mut padded_seq = seq;
    let padded_seq_owned;
    if append_end {
        padded_seq_owned = format!("{}{}{}", "$".repeat(k), seq, "$".repeat(k));
        padded_seq = &padded_seq_owned;
    }
    
    let padded_seq_bytes = padded_seq.as_bytes();

    let mut kmer_vec: Vec<KMer> = Vec::with_capacity(seq.len() + k + 1);
    let len = padded_seq_bytes.len();

    // Init with the first k-1 nucleotides
    let mut rolling_kmer: KMer = KMer::new(k);

    for i in 0..k - 1 {
        // update the k-mers
        rolling_kmer._update_with_nuc(padded_seq_bytes[i]);
    }

    // Iterate through the rest of the string
    for i in k-1..len {
        let nuc_byte = padded_seq_bytes[i];

        // update the k-mers
        rolling_kmer._update_with_nuc(nuc_byte);

        // push the k-mers to the vector
        kmer_vec.push(rolling_kmer.clone());
    }

    kmer_vec
}

/**
 * @brief Convert a quality score to log probability of being correct.
 * [FIXME] Make this more numerically stable.
 * 
 * @param qual The quality score in the fastq format
 * @returns The log probability of being correct
 */
pub fn _qual_to_log_prob(qual: u8) -> f64 {
    let q_score = (qual as i32) - 33;
    let prob_error = 10f64.powf(-(q_score as f64) / 10.0);
    let prob_correct = 1.0 - prob_error;
    prob_correct.ln()
}


/**
 * @brief Convert a quality score to probability of being correct.
 * [FIXME] Make this more numerically stable.
 * 
 * @param qual The quality score in the fastq format
 * @returns The probability of being correct
 */
pub fn qual_to_weight(qual: u8) -> f64 {
    let q_score = (qual as i32) - 33;
    let prob_error = 10f64.powf(-(q_score as f64) / 10.0);
    1.0 - prob_error
}


pub fn qual_edge_kmers(
    qual: &String,
    k: usize
) -> (f64, f64) {
    // If the string is shorter than the k-mer length, return
    if qual.len() < k {
        return (0.0, 0.0);
    }

    let qual_bytes = qual.as_bytes();

    let first_kmer_qual: f64 = (0..k).map(|i| _qual_to_log_prob(qual_bytes[i])).sum::<f64>().exp();
    let last_kmer_qual: f64 = (qual.len() - k..qual.len()).map(|i| _qual_to_log_prob(qual_bytes[i])).sum::<f64>().exp();

    (first_kmer_qual, last_kmer_qual)
}

/**
 * @brief Compute the quality of the k-mers, by summing up the log probabilities
 * of each base being correct.
 * 
 * @param qual The quality scores in the fastq format
 * @param k The length of the k-mers
 * @returns A vector of quality of k-mers
 */
pub fn qual_to_kmer_vec(
    qual: &String,
    k: usize
) -> Vec<f64> {
    // If the string is shorter than the k-mer length, return
    if qual.len() < k {
        return Vec::new();
    }

    let mut kmer_vec: Vec<f64> = Vec::with_capacity(qual.len() - k + 1);

    let qual_vec: Vec<f64> = qual.bytes().map(|q| _qual_to_log_prob(q)).collect();

    // Init with the first k-1 nucleotides
    let mut rolling_kmer: f64 = 0.0;

    for i in 0..k - 1 {
        rolling_kmer += qual_vec[i];
    }

    // Iterate through the rest of the string
    for i in k-1..qual_vec.len() {
        rolling_kmer += qual_vec[i];
        if i >= k {
            rolling_kmer -= qual_vec[i - k];
        }

        // push the k-mers to the vector
        kmer_vec.push(rolling_kmer.exp());
    }

    kmer_vec
}


/**
 * @brief Convert a k-mer string to a u32 index
 */
pub fn kmer_to_index(kmer: &str) -> u32 {
    // use BYTE_TO_SEQ to convert the k-mer to a u32
    let mut res: u32 = 0;
    for b in kmer.bytes() {
        let nuc = BYTE_TO_SEQ[b as usize] as u32;
        res = (res << 2) | nuc;
    }
    res
}

/**
 * @brief Convert a k-mer index to a string
 * 
 * @param kmer The k-mer index
 * @param k The length of the k-mer
 * @returns The k-mer string
 */
pub fn index_to_kmer(kmer: u32, k: usize) -> String {
    // use SEQ_TO_BYTE to convert the k-mer to a string
    let mut kmer_str = String::new();
    for i in (0..k).rev() {
        let nuc = kmer >> (i * 2) & 3;
        kmer_str.push(SEQ_TO_BYTE[nuc as usize] as char);
    }
    kmer_str
}