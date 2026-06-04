use crate::cmdline::*;
use crate::markov_chain::*;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::{mpsc, Arc};
use std::thread;


pub fn consensus(args: Cli) {
    let num_threads = args.threads.unwrap_or_else(rayon::current_num_threads);

    if args.format.to_lowercase() == "microsoft" {
        parse_clusters_file_microsoft(&args, num_threads).unwrap();
    } else if args.format.to_lowercase() == "dna_storage_toolkit" {
        parse_clusters_file_dna_storage_toolkit(&args, num_threads).unwrap();
    } else {
        panic!("Unsupported input format: {}. Currently supports 'microsoft' and 'dna_storage_toolkit'.", args.format);
    }
}

/// A cluster item: either a non-empty cluster (Some) or a consecutive-separator
/// placeholder that should emit an empty line (None).
type ClusterItem = Option<Vec<String>>;

/// Result for a single cluster: (sequence, k, path_weight, confidence).
type ClusterResult = Option<(String, usize, f64, f64)>;

/// Run the parallel pipeline over an iterator of cluster items.
///
/// - A bounded work channel provides backpressure so the reader never loads
///   the whole file into memory.
/// - `num_threads` worker threads each pull from the channel independently.
/// - An ordering buffer in the calling thread re-sequences results and prints them.
fn run_pipeline(
    kmc: Arc<KMarkovChain>,
    clusters: impl Iterator<Item = ClusterItem> + Send + 'static,
    length: usize,
    num_threads: usize,
    output_path: Option<String>,
) {
    // Bound = 2× workers: reader blocks when workers are busy, limiting memory.
    let buffer = num_threads * 2;
    let (work_tx, work_rx) = mpsc::sync_channel::<(usize, ClusterItem)>(buffer);
    let (result_tx, result_rx) = mpsc::channel::<(usize, ClusterResult)>();

    // Shared work queue consumed by all workers.
    let work_rx = Arc::new(std::sync::Mutex::new(work_rx));

    // Spawn worker threads.
    let workers: Vec<_> = (0..num_threads)
        .map(|_| {
            let work_rx = work_rx.clone();
            let result_tx = result_tx.clone();
            let kmc = kmc.clone();
            thread::spawn(move || loop {
                let item = work_rx.lock().unwrap().recv();
                match item {
                    Ok((idx, Some(reads))) => {
                        let (seq, k, weight, confidence) = kmc.find_consensus(length, &reads);
                        result_tx.send((idx, Some((seq, k, weight, confidence)))).unwrap();
                    }
                    Ok((idx, None)) => {
                        result_tx.send((idx, None)).unwrap();
                    }
                    Err(_) => break,
                }
            })
        })
        .collect();

    // Drop our copy of result_tx so the channel closes when all workers finish.
    drop(result_tx);

    // Feed clusters into the work channel from a dedicated reader thread.
    // sync_channel blocks here when the buffer is full → backpressure.
    let reader = thread::spawn(move || {
        for (idx, cluster) in clusters.enumerate() {
            work_tx.send((idx, cluster)).unwrap();
        }
        // work_tx drops here, signalling workers to stop.
    });

    // Collect results and reorder them before printing.
    let mut next_idx: usize = 0;
    let mut out_of_order: HashMap<usize, ClusterResult> = HashMap::new();

    let mut csv_writer = output_path.map(|path| {
        let file = File::create(&path).expect("Failed to create output file");
        let mut w = BufWriter::new(file);
        writeln!(w, "read_id,reconstruction_result,k,path_weight,confidence").unwrap();
        w
    });
    let mut read_id: usize = 1;

    for (idx, result) in result_rx {
        out_of_order.insert(idx, result);
        while let Some(r) = out_of_order.remove(&next_idx) {
            match r {
                Some((ref seq, k, weight, confidence)) => {
                    if let Some(ref mut w) = csv_writer {
                        writeln!(w, "{},{},{},{:.4},{:.6}", read_id, seq, k, weight, confidence).unwrap();
                    } else {
                        println!("{}", seq);
                    }
                    read_id += 1;
                }
                None => {
                    if csv_writer.is_none() {
                        println!();
                    }
                }
            }
            next_idx += 1;
        }
    }

    reader.join().unwrap();
    for w in workers {
        w.join().unwrap();
    }
}


fn parse_clusters_file_dna_storage_toolkit(
    args: &Cli,
    num_threads: usize,
) -> Result<(), Box<dyn std::error::Error>> {

    let kmc = Arc::new(KMarkovChain::new(
        args.k_min as usize, args.k_max as usize,
        args.beam_width as usize, args.alpha as f64,
        !args.single_sided, args.debug,
    ));

    for path in &args.files {
        let file = File::open(path)?;
        let mut lines = BufReader::new(file).lines();

        let num_clusters: usize = lines
            .next()
            .ok_or("Missing number of clusters")??
            .trim()
            .parse()?;

        // Build a lazy iterator that parses one cluster at a time from the file.
        let num_reads_limit = args.num_reads;
        let cluster_iter = (0..num_clusters).map(move |_| -> ClusterItem {
            let n_reads: usize = lines
                .next()
                .expect("Missing number of reads for cluster")
                .expect("IO error")
                .trim()
                .parse()
                .expect("Invalid read count");

            let mut reads = Vec::with_capacity(n_reads);
            let mut qscores = Vec::with_capacity(n_reads);

            for _ in 0..n_reads {
                let line = lines
                    .next()
                    .expect("Unexpected end of file while reading reads")
                    .expect("IO error");
                let line = line.trim().to_string();
                if line.is_empty() {
                    continue;
                }
                let mut parts = line.splitn(2, ',');
                let read = parts.next().expect("Missing read").trim().to_string();
                let q_score = parts.next().expect("Missing q-score").trim().to_string();
                reads.push(read);
                qscores.push(q_score);
            }

            if let Some(limit) = num_reads_limit {
                reads.truncate(limit);
                qscores.truncate(limit);
            }

            Some(reads)
        });

        run_pipeline(kmc.clone(), cluster_iter, args.length, num_threads, args.output_path.clone());
    }

    Ok(())
}


fn parse_clusters_file_microsoft(
    args: &Cli,
    num_threads: usize,
) -> Result<(), Box<dyn std::error::Error>> {

    let kmc = Arc::new(KMarkovChain::new(
        args.k_min as usize, args.k_max as usize,
        args.beam_width as usize, args.alpha as f64,
        !args.single_sided, args.debug,
    ));

    for path in &args.files {
        let file = File::open(path)?;
        let lines = BufReader::new(file).lines();
        let separator = args.separator.clone();

        // Build a lazy iterator that yields one ClusterItem per separator boundary.
        let cluster_iter = MicrosoftClusterIter::new(lines, separator);

        run_pipeline(kmc.clone(), cluster_iter, args.length, num_threads, args.output_path.clone());
    }

    Ok(())
}


/// Iterator that yields one `ClusterItem` per cluster boundary in Microsoft format.
struct MicrosoftClusterIter {
    lines: std::io::Lines<BufReader<File>>,
    separator: String,
    /// Lines buffered from the previous iteration step.
    pending: Option<Vec<String>>,
    /// Whether we've seen the first separator yet.
    first_cluster: bool,
    done: bool,
}

impl MicrosoftClusterIter {
    fn new(lines: std::io::Lines<BufReader<File>>, separator: String) -> Self {
        Self { lines, separator, pending: None, first_cluster: true, done: false }
    }
}

impl Iterator for MicrosoftClusterIter {
    type Item = ClusterItem;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            // Flush the last cluster if any.
            return None;
        }

        let mut current: Vec<String> = self.pending.take().unwrap_or_default();

        loop {
            match self.lines.next() {
                None => {
                    self.done = true;
                    if current.is_empty() {
                        return None;
                    }
                    return Some(Some(current));
                }
                Some(Err(e)) => panic!("IO error: {}", e),
                Some(Ok(line)) => {
                    if line.trim_start().starts_with(self.separator.as_str()) {
                        if self.first_cluster {
                            self.first_cluster = false;
                            // No cluster yet, keep going.
                        } else if current.is_empty() {
                            // Consecutive separator → emit empty line.
                            return Some(None);
                        } else {
                            // End of a cluster.
                            return Some(Some(current));
                        }
                    } else {
                        let t = line.trim().to_string();
                        if !t.is_empty() {
                            current.push(t);
                        }
                        self.first_cluster = false;
                    }
                }
            }
        }
    }
}
