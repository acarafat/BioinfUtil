"""
Gene Present-Absent Matrix

Let's say you have a set of fasta files where each file contains specific gene sequence from different isolates. 
This script can make your life easier to make a gene present-absent matrix.
"""

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
    for filename in listdir(getcwd()):
        if filename.endswith(suffix):
            gene_list.append(filename)
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
    
if __name__ == '__main__':       
    gene_list = enlist_fasta(argv[1], argv[2])
    gene_pa = {}
    
    for gene in gene_list:
        gene_pa[gene] = enlist_isolates(gene)
    
    siMat = pd.DataFrame.from_dict(gene_pa)
