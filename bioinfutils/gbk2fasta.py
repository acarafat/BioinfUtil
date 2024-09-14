import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

'''
A versatile script to convert a GenBank file to fasta. Extract whole contigs, only CDS nucleotide sequences, or translated CDS sequences. This is a Python script to be used in command line/shell/terminal.

Issues and bug-reports are appreciated: https://github.com/acarafat/gbk2fasta

Written by Arafat Rahman
'''

def gbk2fasta(gbk_file):
    '''
    INPUT: A GenBank file containing one or multiple seq-records
    OUTPUT: A list of all seq-records
    '''
    return list(SeqIO.parse(gbk_file, 'genbank'))


def gbk2ffn(gbk_file):
    '''
    INPUT: A GenBank file path
    OUTPUT: A list of nucleotide sequences for CDS features
    '''
    seq_record_list = []

    for gb_record in SeqIO.parse(gbk_file, 'genbank'):
        for gb_feature in gb_record.features:
            if gb_feature.type == 'CDS':
                locus = gb_feature.qualifiers.get('locus_tag', [''])[0]
                product = gb_feature.qualifiers.get('product', [''])[0]
                seq_record = gb_feature.extract(gb_record)
                seq_record.id = locus
                seq_record.name = ''
                seq_record.description = product
                seq_record_list.append(seq_record)

    print(f'Total {len(seq_record_list)} CDS extracted from {gbk_file}')
    return seq_record_list

def gbk2faa(gbk_file):
    '''
    INPUT: A GenBank file path
    OUTPUT: A list of translated amino acid sequences for CDS features
    '''
    aa_seq_list = []

    for gb_record in SeqIO.parse(gbk_file, 'genbank'):
        for gb_feature in gb_record.features:
            if gb_feature.type == 'CDS':
                locus = gb_feature.qualifiers.get('locus_tag', [''])[0]
                product = gb_feature.qualifiers.get('product', [''])[0]
                aa_seq = gb_feature.extract(gb_record.seq).translate(table=11, cds=True)
                aa_seq_record = SeqRecord(aa_seq, id=locus, description=product, name='')
                aa_seq_list.append(aa_seq_record)

    print(f'Total {len(aa_seq_list)} CDS extracted from {gbk_file}')
    return aa_seq_list


def process_file(gbk_file, output_type, output_dir):
    '''
    Processes a single GenBank file and writes the corresponding output
    '''
    base_name = os.path.splitext(os.path.basename(gbk_file))[0]

    if output_type == "fna":
        all_seq = gbk2fasta(gbk_file)
        suffix = '.fna'
    elif output_type == "ffn":
        all_seq = gbk2ffn(gbk_file)
        suffix = '.ffn'
    elif output_type == "faa":
        all_seq = gbk2faa(gbk_file)
        suffix = '.faa'
    else:
        raise ValueError("Invalid output type specified.")

    output_file = os.path.join(output_dir, base_name + suffix)
    SeqIO.write(all_seq, output_file, 'fasta')
    print(f'File saved to {output_file}')


def process_input(input_path, output_type, output_dir):
    '''
    Determines if the input is a single file or directory and processes accordingly
    '''
    if os.path.isfile(input_path):
        process_file(input_path, output_type, output_dir)
    elif os.path.isdir(input_path):
        for file_name in os.listdir(input_path):
            if file_name.endswith(".gbk") or file_name.endswith(".gbff"):
                gbk_file = os.path.join(input_path, file_name)
                process_file(gbk_file, output_type, output_dir)
    else:
        raise ValueError("Invalid input path. Must be a file or directory.")



def main():
    parser = argparse.ArgumentParser(description='Convert GenBank files to fasta format (fna, ffn, faa).')
    parser.add_argument('--input', '-i',  help='Input GenBank file or directory containing GenBank files.')
    parser.add_argument('--output_type', '-t', choices=['fna', 'ffn', 'faa'], help='Output type: fna (whole contigs),
 ffn (CDS nucleotide sequences), faa (translated CDS sequences).')
    parser.add_argument('--output_dir', '-o', help='Directory to save output files.')

    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    process_input(args.input, args.output_type, args.output_dir)


if __name__ == "__main__":
    main()