## Synopsis

[![DOI](https://zenodo.org/badge/175933936.svg)](https://zenodo.org/badge/latestdoi/175933936)

ToyCOATi is a demonstration workflow for using the OpenFST library to reconstruct the maximum-likelihood alignment between two sequences.
You can review the comments in the code to see how it works in more detail.

## Dependencies

 - [OpenFST](http://openfst.org/)
 - [GNU Make](https://www.gnu.org/software/make/)
 - [R](https://www.r-project.org/) and libraries `stringr`, `seqinr`, `Matrix`, and `data.table`.

## Usage

To align the first two sequences in a fasta file, copy that file into the toycoati/fasta directory, then run `make aln/<filename>`.
This will align the sequences and create a file in the `aln` directory.  

### Example

```bash
git clone https://github.com/reedacartwright/toycoati.git toycoati
cd toycoati
cp path/to/example.fasta ./fasta
make aln/example.fasta
```
