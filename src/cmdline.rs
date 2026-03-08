use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about = "Bidirectional beam search for efficient trace reconstruction.", arg_required_else_help = true, disable_help_subcommand = true)]
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

    #[clap(short='S', help_heading = "ALGORITHM", help = "Use single-sided de Bruijn graph. If not set, a double-sided de Bruijn graph will be used.")]
    pub single_sided: bool,

    #[clap(short, default_value_t = 100, help_heading = "INPUT", help = "Length of the target sequences.")]
    pub length: usize,

    #[clap(short, default_value_t = String::from("==="), help_heading = "INPUT", help = "Separator for clusters in Microsoft format. If the line starts with this separator, it marks the boundary/start of a new cluster.")]
    pub separator: String,

    #[clap(short, help_heading = "INPUT", hidden = true, help = "For testing only: number of reads to use in the reconstruction process. If this parameter is not set, all reads in the cluster will be used.")]
    pub num_reads: Option<usize>,

    #[clap(short, help_heading = "OUTPUT", help = "Output file.")]
    pub output_path: Option<String>,

    #[clap(short, help_heading = "DEBUG", help = "Enable debug mode.")]
    pub debug: bool,

    #[clap(long, help_heading = "INPUT", default_value_t = String::from("microsoft"), help = "Format for input clusters, currently supports 'microsoft' (default) and 'dna_storage_toolkit'.")]
    pub format: String,
}
