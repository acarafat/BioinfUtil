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

## Get Fasta Stat `get_fasta_stat.py`
For a set of fasta files in a directory, retrieve number of nucleotide per fasta file.
This is helpful to get estimate of genome size.

Usage: `python get_fasta_stat.py <directory> <fasta_suffix>`

## blastn_miner.py
Used to mine blast output

## splitFastaSeq.py
Split a file containing multiple fasta sequences into individual files containing multiple fasta file
This is a command line file. That means you have to use it in the terminal. 

Use of prefix for the splitted fasta file is optional.

Usage: `python splitFastaSeq.py inputFasta.fasta ~/output/directory/ prefix`
