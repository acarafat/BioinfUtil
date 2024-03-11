# BioinfUtil
Regularly used bioinformatics script for different purposes:
- Make gene present-absent matrix
- Split Fasta file
- GenBank to Fasta conversion
- Update Fasta file description

## gbk2fasta
A varsatile script to convert a GenBank file to fasta. Extract whole contigs, only CDS nucleotide sequences, or translated CDS sequences. This is a Python script to be used in command line/shell/terminal. 

### Requirement
It requires BioPython. 

### Usage
Use the following prompt to convert .gbk file to .fasta:
`python gbk2fasta.py input.gbk fasta_option`


### Conversion options
There are three options: `fna`, `ffn`, and `faa`.
- `fna` extracts all GenBank sequence records.
- `ffn` extracts only CDS records.
- `faa` extracts CDS and translate it into amino acid sequences.
  
### Output
`input.gbk.fasta` file will be created in the same directory of this script file.

## Gene Present-Absent Matrix `gene_pa_matrix.py`
Let's say you have a set of fasta files where each file contains specific gene sequence from different isolates. 
This script can make your life easier to make a gene present-absent matrix. 
It can also count length of the genes in #nucleotide
It will make a csv file containing the gene present-absent matrix.

Options: 
a: just the presence-absence matrix
b: also count gene length for that genome in the presence-absence matrix

Usage in command line: `python gene_pa_matrix.py <option> <~/dir/contianing/fasta> <fasta_suffix>`

## Update fasta description `update_fasta_desc.py`

What it does:
- Use provided pattern to split fasta description
- Keep the expected substring from provided position
- Only keep longest sequence if there are multiple entry for same seqID

Usage:
`python update_fasta_desc.py <~/path/to/input/fasta> <fasta_suffix> <split_symbol> <position_to_keep> <updated_fasta_suffix>`

For example, let's say a fasta description is the following:
`>05LoS16R10_36_00608__05LoS16R10_36`
Here the split symbol will be `__` and position to keep will be 1 (count started from 0). Therefore, it will only keep `05LoS16R10_36`

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

## gbk_reverse_complement.py
The provided Python script is designed to reverse complement a GenBank file, including its DNA sequence and gene features. The script utilizes the Biopython library for bioinformatics tasks. This manual will guide you through the usage of the script and provide step-by-step instructions.

Ensure you have a GenBank file (with a .gb or .gbk extension) that you want to reverse complement.

```python reverse_complement_genbank.py --input input.gb --output output_reversed.gb```

Replace input.gb with the path to your input GenBank file and output_reversed.gb with the desired output file name.

Script Options The script accepts the following command-line arguments:
- --input or -i: Path to the input GenBank file (required).
- --output or -o: Path to the output GenBank file (required).

Ensure that you provide both input and output file paths when running the script.





