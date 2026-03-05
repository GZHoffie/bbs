use std::collections::{HashMap, HashSet};

use crate::seeding::*;
use crate::types::*;

pub struct KMerProfile {
    // Hash map that stores all the k-mer counts
    pub kmer_profile: HashMap<KMer, (f64, f64)>, // that maps a k-mer to the quality of its last base and count
    pub kmer_set: HashSet<KMer>, // set of all kmers that appear in the reads
}


pub struct KMarkovChain {
    // Range of k
    pub k_min: usize,
    pub k_max: usize,

    // Number of paths to keep in the beam search
    pub b: usize,

    // Pseudo-count for the kmer profile
    pub alpha: f64,

    // Debug mode
    pub debug: bool,
}

impl KMarkovChain {
    pub fn new(k_min: usize, k_max: usize, b: usize, alpha: f64, debug: bool) -> Self {
        KMarkovChain { k_min, k_max, b, alpha, debug }
    }


    /**
     * Select k such that k-mers in each read are unique.
     */
    fn select_k(&self, reads_vec: &Vec<String>) -> usize {
        // flag checking if all k-mers in the read are unique
        let mut flag: Vec<bool> = vec![false; reads_vec.len()];
        
        // gradually increase k
        let mut k = self.k_min;
        while k <= self.k_max {
            for i in 0..reads_vec.len() {
                if flag[i] {
                    // skip if the read is already marked as passed
                    continue;
                }
                let kmer_vec_i = seq_to_kmer_vec(&reads_vec[i], k, false);
                // check if the k-mers in the read are unique
                let mut kmer_set: HashMap<KMer, usize> = HashMap::new();
                for kmer in kmer_vec_i.iter() {
                    if kmer_set.get(kmer).is_some() {
                        // if the kmer is already in the set, mark it as false
                        //println!("Read {} has non-unique k-mer {} at k={}", i, index_to_kmer(*kmer, k), k);
                        flag[i] = false;
                        break;
                    } else {
                        // if the kmer is not in the set, add it to the set
                        kmer_set.insert(kmer.clone(), 1);
                        flag[i] = true;
                    }
                }
            }
            if flag.iter().all(|&x| x) {
                // all reads have k-mers of length k
                break;
            }
            k += 1;
        }

        k
    }

    pub fn learn(&self, k: usize, reads_vec: &Vec<String>, quals_vec: Option<&Vec<String>>) -> KMerProfile {

        let mut kmer_profile: HashMap<KMer, (f64, f64)> = HashMap::new();
        let mut kmer_set: HashSet<KMer> = HashSet::new();

        // learn the k-mers from the reads
        if quals_vec.is_none() || reads_vec.len() != quals_vec.unwrap().len() {
            // without quality scores
            for read in reads_vec.iter() {
                let kmer_vec = seq_to_kmer_vec(read, k, true);
                for kmer in kmer_vec.iter() {
                    kmer_set.insert(*kmer);
                    if kmer_profile.get(kmer).is_none() {
                        kmer_profile.insert(*kmer, (1.0, 1.0));
                    } else {
                        kmer_profile.insert(*kmer, (kmer_profile.get(kmer).unwrap().0 + 1.0, kmer_profile.get(kmer).unwrap().1 + 1.0));
                    }
                }
            }
            //for kmer in kmer_profile.keys() {
                // print the content
                //println!("Kmer: {}, Count: {:?}", kmer.to_string(), kmer_profile.get(kmer).unwrap());
            //}
            return KMerProfile {
                kmer_profile,
                kmer_set,
            }
        }
        let quals_vec_some = quals_vec.unwrap();
        for (read, qual) in reads_vec.iter().zip(quals_vec_some.iter()) {
            //println!("Read: {}, Qual: {}", read, qual);
            // length = l - k + 1 + 2k = l + k + 1
            let kmer_vec = seq_to_kmer_vec(read, k, true);

            // for the other k-mers, use the quality scores
            let qual_bytes = qual.as_bytes();

            let qual_weights: Vec<f64> = qual_bytes.iter().map(|&q| qual_to_weight(q)).collect();
            // append k 1.0s to both ends of qual weights to represent the padded '$' characters.
            // length = l + 2k
            //println!("Qual len: {}, Read len: {}", qual_weights.len(), read.len());
            assert!(qual.len() == read.len());
            let mut qual_weights_padded: Vec<f64> = vec![1.0; k];
            qual_weights_padded.extend(qual_weights);
            qual_weights_padded.extend(vec![1.0; k]);

            assert!(qual_weights_padded.len() == read.len() + 2 * k);

            //println!("Padded qual weights: {:?}", qual_weights_padded);

            
            for i in 0..kmer_vec.len() {
                let kmer = &kmer_vec[i];
                kmer_set.insert(*kmer);

                // update kmer_profile with the quality of the first and last base of the k-mer
                let last_base_weight = qual_weights_padded[i + k - 1];
                let first_base_weight = qual_weights_padded[i];

                if kmer_profile.get(kmer).is_none() {
                    kmer_profile.insert(*kmer, (first_base_weight, last_base_weight));
                } else {
                    kmer_profile.insert(*kmer, (kmer_profile.get(kmer).unwrap().0 + first_base_weight, kmer_profile.get(kmer).unwrap().1 + last_base_weight));
                }


                // update the weight of the neighbors of the k-mer 
                if last_base_weight < 1.0 {
                    for neighbor in kmer.neighbor_kmers_last_base() {
                        if neighbor == *kmer {
                            continue;
                        }
                        if kmer_profile.get(&neighbor).is_none() {
                            kmer_profile.insert(neighbor, (0.0, (1.0 - last_base_weight) / 3.0));
                        } else { 
                            kmer_profile.insert(neighbor, (kmer_profile.get(&neighbor).unwrap().0, kmer_profile.get(&neighbor).unwrap().1 + (1.0 - last_base_weight) / 3.0));
                        }
                    }
                }

                if first_base_weight < 1.0 {
                    for neighbor in kmer.neighbor_kmers_first_base() {
                        if neighbor == *kmer {
                            continue;
                        }
                        if kmer_profile.get(&neighbor).is_none() {
                            kmer_profile.insert(neighbor, ((1.0 - first_base_weight) / 3.0, 0.0));
                        } else { 
                            kmer_profile.insert(neighbor, (kmer_profile.get(&neighbor).unwrap().0 + (1.0 - first_base_weight) / 3.0, kmer_profile.get(&neighbor).unwrap().1));
                        }
                    }
                }
            } 
        }

        //println!("{:?}", kmer_profile);
        //for (kmer, count) in kmer_profile.iter() {
        //    println!("Kmer: {}, Count: {}", index_to_kmer(*kmer, k), count);
        //}

        //for kmer in kmer_profile.keys() {
            // print the content
            //println!("Kmer: {}, Count: {:?}", kmer.to_string(), kmer_profile.get(kmer).unwrap());
        //}
        KMerProfile {
            kmer_profile,
            kmer_set,
        }
    }

    
    /* BEAM SEARCH UTILITIES */
    fn beam_search_step(
        &self,
        k: usize,
        frontier: &Vec<(Vec<KMer>, f64)>,
        profile: &KMerProfile,
        forward: bool,
        next_base_end: bool,
        ending_kmer: KMer,
    ) -> Vec<(Vec<KMer>, f64)> {
        let mut new_frontier: Vec<(Vec<KMer>, f64)> = Vec::new();

        for (kmer_vec, score) in frontier.iter() {
            let current_kmer = *kmer_vec.last().unwrap();
            let mut total_count: f64 = 0.0;

            // first iteration to find the total count of next kmers
            // All the k-mers that can be the next k-mer in the de Bruijn graph
            let next_kmers = if forward {
                if next_base_end {
                    current_kmer.next_kmers_at_end(ending_kmer)
                } else {
                    current_kmer.next_kmers()
                }
            } else {
                if next_base_end {
                    current_kmer.prev_kmers_at_end(ending_kmer)
                } else {
                    current_kmer.prev_kmers()
                }
            };

            // All the k-mers that can be the next k-mer in the consensus
            // [TODO] Take ECC into consideration
            let next_valid_kmers = if next_base_end {
                // only consider the next kmers that has the ending marker
                if forward {
                    vec![current_kmer.next_kmer_end()]
                } else {
                    vec![current_kmer.prev_kmer_end()]
                }
            } else {
                next_kmers.clone()
            };

            //for next_kmer in &next_kmers {
            //    println!("Next kmer: {}, profile: {:?}", next_kmer.to_string(), profile.kmer_profile.get(next_kmer));
            //}

            for next_kmer in &next_kmers {
                if profile.kmer_set.get(next_kmer).is_none() {
                    continue;
                }
                if profile.kmer_profile.get(next_kmer).is_some() {
                    total_count += if forward { profile.kmer_profile.get(next_kmer).unwrap().1 } else { profile.kmer_profile.get(next_kmer).unwrap().0 };

                }
            }
            
            //println!("Total count: {}", total_count);



            // second iteration to find the new frontier
            for next_kmer in &next_valid_kmers {
                if profile.kmer_set.get(next_kmer).is_none() {
                    continue;
                }
                let count = if profile.kmer_profile.get(next_kmer).is_some() {
                    if forward { profile.kmer_profile.get(next_kmer).unwrap().1 } else { profile.kmer_profile.get(next_kmer).unwrap().0 }
                } else {
                    0.0
                };
                if count > 0.0 {
                    let new_score = score + ((count + 1.0 * self.alpha) / (total_count + 4.0 * self.alpha)).log2();
                    
                    // append the new kmer to the kmer_vec
                    let mut kmer_vec_new = kmer_vec.clone();
                    kmer_vec_new.push(*next_kmer);
                    new_frontier.push((kmer_vec_new, new_score));
                }
            }

            // keep the top b paths in the new frontier
            // [TODO] Use O(n) partitioning algorithm to find the top b paths
            new_frontier.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            new_frontier.truncate(self.b);
        }

        //println!("New frontier: {:?}", new_frontier);
        //for (kmer_vec, score) in new_frontier.iter() {
        //    println!("  Kmer vec: {}, score: {}", self.path_to_sequence(kmer_vec, k, forward), score);
        //}

        new_frontier
    }


    pub fn beam_search(
        &self,
        k: usize,
        profile: &KMerProfile,
        target_length: usize,
        forward: bool
    ) -> (Vec<KMer>, f64) {
        // set the initial frontier as <'$' * k, 0.0>
        let init_kmer = KMer::from_str(k, &"$".repeat(k));
        let mut frontier: Vec<(Vec<KMer>, f64)> = vec![(vec![init_kmer], 0.0)];

        if self.debug {
            
            for (kmer_vec, score) in frontier.iter() {
                println!("  Kmer vec: {}, score: {}", self.path_to_sequence(kmer_vec, k, forward), score);
            }
        }

        let goal_kmer = KMer::from_str(k, &"$".repeat(k));

        // forward search
        for step in 1..=target_length+k {
            if self.debug {
                println!("Step {}: ", step);
            }
            let new_frontier = self.beam_search_step(k, &frontier, profile, forward, step > target_length, goal_kmer);
            if self.debug {
                for (kmer_vec, score) in new_frontier.iter() {
                
                    println!("  Kmer vec: {}, score: {}", self.path_to_sequence(kmer_vec, k, forward), score);
                }
            }
            if new_frontier.is_empty() {
                // no more paths to explore, failure.
                break;
            }
            frontier = new_frontier;
        }
        
        // find the best path that ends with end_kmer
        let mut best_path = Vec::new();
        let mut best_score = f64::MIN;
        let mut found_target = false;
        for (kmer_vec, score) in frontier.iter() {
            let match_target = (kmer_vec.last().unwrap() == &goal_kmer);
            if ((!found_target) && match_target) || (((found_target && match_target) || (!found_target && !match_target)) && (*score as f64 > best_score)) {
                best_path = kmer_vec.clone();
                best_score = *score as f64;
                if match_target {
                    found_target = true;
                }
            }
        }

        if self.debug {
            println!("Best path: {}, score: {}", self.path_to_sequence(&best_path, k, forward), best_score);
        }

        // return the best path
        (best_path, best_score)
    }

    pub fn path_to_sequence(&self, kmer_vec: &Vec<KMer>, k: usize, forward: bool) -> String {
        if kmer_vec.is_empty() {
            return String::new();
        }

        if forward {
            let mut sequence: String = kmer_vec[0].to_string();
    
            for i in 1..kmer_vec.len() {
                let base = kmer_vec[i].last_nuc_as_u32();
                sequence.push(SEQ_TO_BYTE[base as usize] as char);
            }
            sequence
        } else {
            let mut sequence: String = kmer_vec[0].to_string();
    
            for i in 1..kmer_vec.len() {
                let base = kmer_vec[i].first_nuc_as_u32();
                sequence.insert(0, SEQ_TO_BYTE[base as usize] as char);
            }
            sequence
        }
    }

    // function that uses all the above to find the consensus sequence
    pub fn find_consensus(&self, target_length: usize, reads_vec: &Vec<String>, quals_vec: Option<&Vec<String>>) -> (String, f64) {
        // find a suitable k
        let k = self.select_k(&reads_vec);
        let chosen_k = k + 1;
        let max_steps = target_length - chosen_k + 1;

        // learn the kmer profile
        let profile = self.learn(chosen_k, reads_vec, quals_vec);

        // perform beam search in both direction
        let (best_path_forward, best_score_forward) = self.beam_search(chosen_k, &profile, target_length, true);

        let (best_path_backward, best_score_backward) = self.beam_search(chosen_k, &profile, target_length, false);

        
        let (best_path, best_score) = if (best_score_forward >= best_score_backward && best_path_forward.len() == best_path_backward.len()) ||(best_path_forward.len() > best_path_backward.len()) {
            (self.path_to_sequence(&best_path_forward, chosen_k, true), best_score_forward)
        } else {
            (self.path_to_sequence(&best_path_backward, chosen_k, false), best_score_backward)
        };

        (best_path, best_score)
    }
}