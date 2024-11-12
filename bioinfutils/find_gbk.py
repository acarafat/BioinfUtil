#!/usr/bin/env python

import argparse
from Bio import SeqIO

def load_genes(file):
    with open(file, "r") as f:
        genes = [line.strip() for line in f.readlines()]
    return genes

def check_genes_in_genbank(gene_list, genbank_file):
    found_genes = []
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene" and "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
                if gene_name in gene_list:
                    print(f"{genbank_file}, {gene_name}")
                    found_genes.append(gene_name)
    return found_genes

def main(args=None):

    parser = argparse.ArgumentParser(description="Check if genes are present in a GenBank file")
    
    parser.add_argument("-g", "--genbank", required=True, help="Input GenBank file")
    parser.add_argument("-l", "--gene_list", required=True, help="Text file with list of gene names (one per line)")
    
    args = parser.parse_args(args)

    gene_list = load_genes(args.gene_list)
    check_genes_in_genbank(gene_list, args.genbank)

if __name__ == "__main__":
    main()