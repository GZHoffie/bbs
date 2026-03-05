pub mod seeding;
pub mod types;
pub mod consensus;
pub mod markov_chain;
pub mod cmdline;

/*
Unit tests

*/

#[cfg(test)]
mod tests {
    use needletail::kmer;

    use super::*;
    use crate::seeding::*;

    #[test]
    fn test_kmer_from_str() {
        let kmer_str = "ACGT";
        let kmer = KMer::from_str(4, kmer_str);
        assert_eq!(kmer.to_string(), "ACGT");
    }

    #[test]
    fn test_next_kmer() {
        let kmer_str = "$$$$";
        let kmer = KMer::from_str(4, kmer_str);
        let next_kmer = kmer.next_kmers();
        println!("Next kmers of {}: {:?}", kmer_str, next_kmer.iter().map(|k| k.to_string()).collect::<Vec<String>>());
        assert_eq!(1, 2);

    }
}