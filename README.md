# BioinfUtil
Regularly used bioinformatics script for different purposes:
- Make gene present-absent matrix
- Split Fasta

## gene_pa_matrix.py

Let's say you have a set of fasta files where each file contains specific gene sequence from different isolates. 
This script can make your life easier to make a gene present-absent matrix.

### Command line usages:
`python3 ~/path/to/fasta/files fasta`

## blastn_miner.py
Used to mine blast output

## splitFastaSeq.py
Split a file containing multiple fasta sequences into individual files containing multiple fasta file

## Usage
This is a command line file. That means you have to use it in the terminal. You can use the following prompt:
`python splitFastaSeq.py inputFasta.fasta ~/output/directory/ prefix`

Use of prefix for the splitted fasta file is optional.

### Requirement
This requires Biopython
