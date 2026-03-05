use clap::Parser;


use bbs::consensus;
use bbs::cmdline::*;


fn main() {
    let cli = Cli::parse();
    match cli.mode {
        Mode::Consensus(consensus_args) => consensus::consensus(consensus_args),    
    }
}

