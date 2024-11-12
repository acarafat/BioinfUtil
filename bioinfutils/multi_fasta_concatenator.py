from Bio import SeqIO
from collections import defaultdict

import os
import argparse

'''
Usage: python3 multi_fasta_concatenator.py input_directory output_file.fasta

Description:
Concatenate sequences from multiple FASTA files by strain. The script takes an input directory containing FASTA files, where each file represents an individual strain with multiple genes extracted from the strain's genome. It concatenates the gene sequences for each strain, filling gaps for missing genes, and writes the concatenated sequences into a single FASTA file.

Arguments:
  input_directory    Path to the input directory containing FASTA files.
  output_file.fasta  Path to the output FASTA file where concatenated sequences will be written.

Example:
  python script_name.py input_directory output_file.fasta

Note:
  - The input directory should contain FASTA files, each representing an individual strain.
  - The output file will contain concatenated sequences for each strain, where gene order is maintained across all strains.
  - If a strain does not have a particular gene, gaps will be used to maintain alignment.
'''

def parse_fasta_files(input_dir):
    """
    Parse multiple FASTA files from an input directory and store sequences by strain name and gene.

    Parameters:
        input_dir (str): Path to the input directory containing FASTA files.

    Returns:
        dict: Dictionary containing sequences grouped by strain name and gene.
    """
    sequences_by_strain_gene = defaultdict(dict)

    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fasta') or file_name.endswith('.fa'):
            fasta_file = os.path.join(input_dir, file_name)
            strain_name = file_name.split('.')[0]  # Assuming file names are strain names

            for record in SeqIO.parse(fasta_file, 'fasta'):
                gene_name = record.id
                gene_sequence = str(record.seq)
                sequences_by_strain_gene[strain_name][gene_name] = gene_sequence


    return sequences_by_strain_gene

def concatenate_sequences(sequences_by_strain_gene):
    """
    Concatenate sequences by strain, filling gaps for missing genes.
    
    Parameters:
        sequences_by_strain_gene (dict): Dictionary containing sequences grouped by strain name and gene.
    
    Returns:
        dict: Dictionary containing concatenated sequences by strain.
    """
    concatenated_sequences = {}
    all_genes = set()
    
    # Collect all genes present in any strain
    for strain_sequences in sequences_by_strain_gene.values():
        all_genes.update(strain_sequences.keys())
    
    # Concatenate sequences for each strain
    for strain_name, strain_sequences in sequences_by_strain_gene.items():
        concatenated_sequence = ''
        for gene in all_genes:
            if gene in strain_sequences:
                concatenated_sequence += strain_sequences[gene]
            else:
                concatenated_sequence += '-' * len(next(iter(strain_sequences.values())))
        concatenated_sequences[strain_name] = concatenated_sequence
    
    return concatenated_sequences

def write_concatenated_sequences(concatenated_sequences, output_file):
    """
    Write concatenated sequences to a FASTA file.
    
    Parameters:
        concatenated_sequences (dict): Dictionary containing concatenated sequences by strain.
        output_file (str): Path to the output FASTA file.
    """
    with open(output_file, 'w') as f:
        for strain_name, sequence in concatenated_sequences.items():
            f.write(f'>{strain_name}\n{sequence}\n')

def main(args=None):
    parser = argparse.ArgumentParser(description='Concatenate sequences from multiple FASTA files by strain')
    parser.add_argument('--input', '-i', help='Input directory containing FASTA files')
    parser.add_argument('--output', '-o', help='Output FASTA file name')
    args = parser.parse_args(args)

    sequences_by_strain_gene = parse_fasta_files(args.input_dir)
    concatenated_sequences = concatenate_sequences(sequences_by_strain_gene)
    write_concatenated_sequences(concatenated_sequences, args.output_file)

if __name__ == "__main__":
    main()
