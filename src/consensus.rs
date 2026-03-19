

use crate::cmdline::*;
use crate::markov_chain::*;

use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};


pub fn consensus(args: ConsensusArgs) {
    // [TODO] Allow other format of input
    if args.format.to_lowercase() == "microsoft" {
        parse_clusters_file_microsoft(&args).unwrap();
        return;
    } else if args.format.to_lowercase() == "dna_storage_toolkit" {
        parse_clusters_file_dna_storage_toolkit(&args).unwrap();
    } else {
        panic!("Unsupported input format: {}. Currently supports 'microsoft' and 'dna_storage_toolkit'.", args.format);
    }
}

fn parse_clusters_file_dna_storage_toolkit(
    args: &ConsensusArgs,
) -> Result<(), Box<dyn std::error::Error>> {

    let kmc = KMarkovChain::new(args.k_min as usize, args.k_max as usize, args.beam_width as usize, args.alpha as f64, !args.single_sided, args.debug);


    for path in &args.files {
        let file = File::open(path)?;
        let mut lines = BufReader::new(file).lines();

        // First line: number of clusters
        let num_clusters: usize = lines
            .next()
            .ok_or("Missing number of clusters")??
            .trim()
            .parse()?;

        //println!("Number of clusters: {}", num_clusters);

        //let mut all_reads: Vec<Vec<String>> = Vec::with_capacity(num_clusters);
        //let mut all_qscores: Vec<Vec<f64>> = Vec::with_capacity(num_clusters);

        for cluster_idx in 0..num_clusters {
            // Line with number of reads in this cluster
            let n_reads_line = lines
                .next()
                .ok_or("Missing number of reads for cluster")??;
            let n_reads: usize = n_reads_line.trim().parse()?;

            let mut reads = Vec::with_capacity(n_reads);
            let mut qscores = Vec::with_capacity(n_reads);

            for _ in 0..n_reads {
                let line = lines
                    .next()
                    .ok_or("Unexpected end of file while reading reads")??;
                let line = line.trim();
                if line.is_empty() {
                    continue; // or return Err(...) if empty lines are invalid
                }

                // Expect: "ReadID,qscore"
                let mut parts = line.splitn(2, ',');
                let read = parts
                    .next()
                    .ok_or("Missing read")?
                    .trim()
                    .to_string();

                let q_score: String;
                q_score = parts
                    .next()
                    .ok_or("Missing q-score")?
                    .trim()
                    .to_string();
                

                reads.push(read);
                qscores.push(q_score);
            }

            if let Some(num_reads) = args.num_reads {
                reads.truncate(num_reads);
                qscores.truncate(num_reads);
            }

            
            
            //println!("Cluster {}: {} reads", cluster_idx + 1, n_reads);
            //println!("Reads: {:?}", reads);
            //println!("Q-scores: {:?}", qscores);
            // find consensus for this cluster
            let (consensus_sequence, score) = kmc.find_consensus(args.length, &reads);

            // delete the sequences of '$' at the beginning and end of the consensus sequence
            let consensus_sequence = consensus_sequence.trim_matches('$').to_string();
            println!("{}", consensus_sequence);
        }
    }

    Ok(())

}


fn parse_clusters_file_microsoft(
    args: &ConsensusArgs,
) -> Result<(), Box<dyn std::error::Error>> {

    let kmc = KMarkovChain::new(args.k_min as usize, args.k_max as usize, args.beam_width as usize, args.alpha as f64, !args.single_sided, args.debug);


    for path in &args.files {
        let file = File::open(path)?;
        let mut lines = BufReader::new(file).lines();



        // Split into clusters: a line starting with "===" marks a boundary/start of a new cluster.
        // Collect non-separator, non-empty lines between separators as one cluster.
        
        let mut first_cluster = true;
        let mut current_cluster: Vec<String> = Vec::new();
        //let mut qscores: Vec<String> = Vec::new();

        //println!("Processing file: {}", path);
        while let Some(line) = lines.next().transpose()? {
            if line.trim_start().starts_with(args.separator.as_str()) {
                if current_cluster.is_empty() {
                    if first_cluster {
                        // Just starting, no cluster to process yet
                        first_cluster = false;
                        continue;
                    } else {
                        // Consecutive separators, skip
                        println!("");
                        continue;
                    }
                } else {
                    let (consensus_sequence, score) = kmc.find_consensus(args.length, &current_cluster);
                    let consensus_sequence = consensus_sequence.trim_matches('$').to_string();
                    println!("{}", consensus_sequence);
                    current_cluster.clear();
                }
            } else {
                let t = line.trim();
                if !t.is_empty() {
                    current_cluster.push(t.to_string());
                }
                first_cluster = false;
            }
        }
        if !current_cluster.is_empty() {
            let (consensus_sequence, score) = kmc.find_consensus(args.length, &current_cluster);
            let consensus_sequence = consensus_sequence.trim_matches('$').to_string();
            println!("{}", consensus_sequence);
        }
    }
        
        
    

    Ok(())

}
