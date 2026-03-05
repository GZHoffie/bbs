use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about = "b3s", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// Find the consensus sequence from the reads
    #[clap(display_order = 1)]
    Consensus(ConsensusArgs),
}

#[derive(Args, Default)]
pub struct ConsensusArgs {
    #[clap(multiple=true, help_heading = "INPUT", help = "Clustered read files.")]
    pub files: Vec<String>,

    #[clap(short='k', default_value_t = 4, help_heading = "ALGORITHM", help ="Minimum value of k for de Bruijn graph construction.")]
    pub k_min: u8,

    #[clap(short='K', default_value_t = 31, help_heading = "ALGORITHM", help ="Maximum value of k for de Bruijn graph construction.")]
    pub k_max: u8,

    #[clap(short, default_value_t = 1.0, help_heading = "ALGORITHM", help ="Value of alpha for smoothing.")]
    pub alpha: f64,

    #[clap(short, default_value_t = 20, help_heading = "ALGORITHM", help ="Beam width used in beam search.")]
    pub beam_width: u8,

    #[clap(short, help_heading = "INPUT", help = "Include quality scores.")]
    pub with_qscores: bool,

    #[clap(short, default_value_t = 100, help_heading = "INPUT", help = "Length of the target sequences.")]
    pub length: usize,

    #[clap(short, help_heading = "INPUT", help = "For testing only: number of reads to use in the reconstruction process. If this parameter is not set, all reads in the cluster will be used.")]
    pub num_reads: Option<usize>,

    #[clap(short, default_value_t = String::new(), help_heading = "OUTPUT", help = "Output file.")]
    pub output_path: String,

    #[clap(short, help_heading = "DEBUG", help = "Enable debug mode.")]
    pub debug: bool,

    #[clap(long, help_heading = "INPUT", help = "Microsoft format for input clusters.")]
    pub microsoft_format: bool,
}
