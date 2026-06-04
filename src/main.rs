use clap::Parser;

use bbs::consensus;
use bbs::cmdline::Cli;

fn main() {
    let cli = Cli::parse();
    consensus::consensus(cli);
}

