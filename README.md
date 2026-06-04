
<p align="center">
  <img src="./logo.png" alt="Logo" width="400"/>
</p>

# Efficient Trace Reconstruction for DNA data storage systems using Bidirectional Beam Search


## Installation

### Download executable

Simply download the executable from the latest release, via the following

```bash
wget https://github.com/GZHoffie/bbs/releases/download/v0.2.0/bbs
```

### Build from source

Alternatively, build skiver from the source code. Install [rust](https://rust-lang.org/tools/install/), and build using

```
https://github.com/GZHoffie/bbs.git
cd bbs

# If default rust install directory is ~/.cargo
cargo install --path . --root ~/.cargo
```

BBS is available by running 

```
bbs
```

## Quick start

Currently, BBS supports 2 input formats,

- `dna_storage_toolkit` format.
- `microsoft` format: similar to the [CNR dataset](https://github.com/microsoft/clustered-nanopore-reads-dataset). Below is a quick demo.


```bash
# Download the Microsoft CNR dataset
git clone https://github.com/microsoft/clustered-nanopore-reads-dataset.git

# Assume that we are still in the build/ directory
# To output the reconstructed sequence directly, use `./bbs -i <input_clusters> -l <read_length> > <output_file_name>`
bbs clustered-nanopore-reads-dataset/Clusters.txt -l 110 > output.txt
```

For more detailed output, use the following command.

```bash
bbs clustered-nanopore-reads-dataset/Clusters.txt -l 110 -o output_verbose.csv
```

In the output `csv` file, there will be 5 fields, indicating the index of the cluster, reconstructed sequence, the value of `k` used, the total path weight, and the confidence value.

```
read_id,reconstruction_result,k,path_weight,confidence
1,ACCATAATGCGTGGGGCCGACCTCGGAATGCGGTCTCCATGCGCGTTTCCTCCAACCTAAGGTAGCCTGTAGTTCATTGGACCTCTGATGGCGCTTATAGAAACCGGGAA,11,-14.9066,0.909951
2,TCGAAGCAGTAGGGCCTACCAAATAGGTTGGTCCTCCGTTGTATCTAAGGATTGAGTTTACCTGGCTTACACGGCAGGTACCGCCAATCTCGTCCGGCTCCGCGGCATCC,8,-32.2539,0.950223
3,AGTTAACGTCCCACGGCGAGGCACTCTTGATCCCCACCTTCAAGAGGTGTACCGGATCATGGAGAACAAGCATACGTCGCACGCACACCATTGGACGGCGAGTGCCGAGT,10,-44.4446,0.853414
```

Use the following for a detailed guide on other input parameters.

```bash
bbs -h
```

## Citation

Gu Z, Xin H, Sharma P, Goh GY, Wong L, Nagarajan N. [Efficient trace reconstruction in DNA storage systems using bidirectional beam search](https://www.cell.com/iscience/fulltext/S2589-0042(25)02052-8). *iScience*. 2025 Oct 21;28(11):113791. doi: 10.1016/j.isci.2025.113791. PMID: 41280695; PMCID: PMC12630026.
