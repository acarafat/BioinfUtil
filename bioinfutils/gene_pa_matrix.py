import argparse
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



def main(args=None):
    parser = argparse.ArgumentParser(description='Convert GenBank files to fasta format (fna, ffn, faa).')
    parser.add_argument('--input', '-i',  help='A directory path that contains multiple fasta files. Each fasta should contains isolates for a specific genes.')
    parser.add_argument('--suffix', '-s', help='Suffix of fasta files')
    parser.add_argument('--option', '-o',  help="Option 'a' to generate presence matrix or 'b' to generate size matrix.")
    args = parser.parse_args(args)

    gene_list = enlist_fasta(args.input, args.suffix)
    
    gene_pa = {}
    if args.option == 'a':
        for gene in gene_list:
            gene_pa[gene] = enlist_entry(gene)
    elif args.option == 'b':
        for gene in gene_list:
            gene_pa[gene] = enlist_entry_size(gene)
    
    siMat = pd.DataFrame.from_dict(gene_pa)
    siMat.to_csv('pa_mat.csv', sep=',')

if __name__ == '__main__':
    main()