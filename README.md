# BioinfUtil
Regularly used bioinformatics script for different purposes:
- Make gene present-absent matrix
- Split Fasta

## Gene Present-Absent Matrix `gene_pa_matrix.py`
Let's say you have a set of fasta files where each file contains specific gene sequence from different isolates. 
This script can make your life easier to make a gene present-absent matrix. 
It can also count length of the genes in #nucleotide
It will make a csv file containing the gene present-absent matrix.

Options: 
a: just the presence-absence matrix
b: also count gene length for that genome in the presence-absence matrix

Usage in command line: `python gene_pa_matrix.py <option> <~/dir/contianing/fasta> <fasta_suffix>`

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
