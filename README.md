
<p align="center">
  <img src="./logo.png" alt="Logo" width="400"/>
</p>

# Efficient Trace Reconstruction for DNA data storage systems using Bidirectional Beam Search


## Installation

### Download executable

Simply download the executable from the latest release, via the following

```bash
wget https://github.com/GZHoffie/bbs/releases/download/v0.0.1/bbs
```

### Build from source

Alternatively, build BBS from the source code. This tool is built with [Seqan3](https://docs.seqan.de/seqan/3-master-user/index.html). To properly build the package, you need to have GCC >= 11.3, G++ and CMake installed.

```bash
git clone https://github.com/GZHoffie/bbs.git
cd bbs
mkdir build
cd build
```

Build the project with the command

```
cmake ../
cmake --build . --target bbs
```

This is going to download all the dependencies, and the executable is then available in `build/bbs`.

## Quick start

Currently, BBS supports input format that is similar to the [CNR dataset](https://github.com/microsoft/clustered-nanopore-reads-dataset). Below is a quick demo.

```bash
# Download the Microsoft CNR dataset
git clone https://github.com/microsoft/clustered-nanopore-reads-dataset.git

# Assume that we are still in the build/ directory
# To output the reconstructed sequence directly, use `./bbs -i <input_clusters> -l <read_length> -s <cluster_separator> > <output_file_name>`
./bbs -i clustered-nanopore-reads-dataset/Clusters.txt -l 110 -s "====" > output.txt
```

The separator is an indicator string that separates the clusters. If we read a line that contains the specified separator, we regard it as a new cluster.

For more detailed output, use the following command.

```bash
./bbs -i clustered-nanopore-reads-dataset/Clusters.txt -l 110 -s "====" -o output.csv > /dev/null
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
./bbs -h
```

## Citation

[[1](https://www.biorxiv.org/content/10.1101/2025.04.16.644694v1)] Gu, Zhenhao, et al. "Efficient trace reconstruction in DNA storage systems using Bidirectional Beam Search." bioRxiv (2025): 2025-04.
