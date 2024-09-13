'''
Gene Present-Absent Matrix

Let's say you have a set of fasta files where each file contains specific gene sequence from different isolates. 
This script can make your life easier to make a gene present-absent matrix. 
It can also count length of the genes in #nucleotide
It will make a csv file containing the gene present-absent matrix.

Options: 
a: just the presence-absence matrix
b: also count gene length for that genome in the presence-absence matrix

Usage in command line: python gene_pa_matrix.py <option> <~/dir/contianing/fasta> <fasta_suffix>
'''

from sys import argv
from Bio import SeqIO
from os import listdir, getcwd
import pandas as pd


def enlist_fasta(dir_fasta, suffix):
    '''
    dir_fasta: A directory path that contains multiple fasta files. Each fasta should contains isolates for a specific genes.
    SUFFIX: Extension of the set of fasta file.
    OUTPUT: List containing name of the gene fasta files in that directory.
    '''
    gene_list = []
    if not dir_fasta.endswith('/'):
        dir_fasta = dir_fasta + '/'
    for filename in listdir(dir_fasta):
        if filename.endswith(suffix):
            gene_list.append(dir_fasta+filename)
    return gene_list

def enlist_entry(gene_fasta):
    '''
    INPUT: fasta file containing gene sequences form different isolates
    OUTPUT: List of isolates that has that gene
    '''
    isolate_list = {}
    for seq_record in SeqIO.parse(gene_fasta, 'fasta'):
        seq_id = seq_record.id
        isolate_list[seq_id] = '+'
    return(isolate_list)

def enlist_entry_size(gene_fasta):
    '''
    INPUT: fasta file containing gene sequences form different isolates
    OUTPUT: List of gene length for isolates that has that gene
    '''
    isolate_list = {}
    for seq_record in SeqIO.parse(gene_fasta, 'fasta'):
        seq_id = seq_record.id
        isolate_list[seq_id] = len(str(seq_record.seq))
    return(isolate_list)



def main():
    gene_list = enlist_fasta(argv[2], argv[3])
    gene_pa = {}
    if argv[1] == 'a':
        for gene in gene_list:
            gene_pa[gene] = enlist_entry(gene)
    elif argv[1] == 'b':
        for gene in gene_list:
            gene_pa[gene] = enlist_entry_size(gene)
    
    siMat = pd.DataFrame.from_dict(gene_pa)
    siMat.to_csv('pa_mat.csv', sep=',')

if __name__ == '__main__':
    main()