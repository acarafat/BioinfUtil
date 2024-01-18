import sys, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

'''
A varsatile script to convert a GenBank file to fasta. Extract whole contigs, only CDS nucleotide sequences, or translated CDS sequences. This is a Python script to be used in command line/shell/terminal.

Issues and bug-reports are appreciated: https://github.com/acarafat/gbk2fasta

Written by Arafat Rahman
'''

def gbk2fasta(gbkFile):
    '''
    INPUT: A GenBank file containing one or multiple seq-record
    OUTPUT: It writes a single fasta file containing all seq-records
    '''
    allSeq = []
    for seq in SeqIO.parse(gbkFile, 'genbank'):
        allSeq.append(seq)

    return allSeq


def gbk2ffn(gbkPath):
    '''
    Iterate through GenBank features and look for CDS, update seq id and seq description
    INPUT example: database/gbk/inoc131.gbk
    '''

    seq_record_list = []

    for gb_record in SeqIO.parse(gbkPath, 'genbank'):
        for gb_feature in gb_record.features:
            if gb_feature.type == 'CDS':
                locus = gb_feature.qualifiers['locus_tag']
                product = gb_feature.qualifiers['product']
                seq_record = gb_feature.extract(gb_record)
                seq_record.id = locus[0]
                seq_record.name = ''
                seq_record.description = product[0]
                seq_record_list.append(seq_record)
    print(f'Total {len(seq_record_list)} CDS extracted')

    return seq_record_list


def gbk2faa(gbkPath):
    '''Iterate through GenBank feature and look for CDS, extract amino acid sequence
    INPUT example: 'database/gbk/inoc131.gbk'
    '''
    aa_seq_list = []

    for gb_record in SeqIO.parse(gbkPath, 'genbank'):
        for gb_feature in gb_record.features:
            if gb_feature.type == 'CDS':
                locus = gb_feature.qualifiers['locus_tag'][0]
                product = gb_feature.qualifiers['product'][0]
                aa_seq = gb_feature.extract(gb_record.seq).translate(table=11, cds=True)
                aa_seqRecord = SeqRecord(aa_seq, id=locus, description=product, name='')
                aa_seq_list.append(aa_seqRecord)
    print(f'Total {len(aa_seq_list)} CDS extracted')

    return aa_seq_list


if __name__ == "__main__":
    inputFile = sys.argv[1]
    outputType = sys.argv[2]
    outputFile = os.path.basename(inputFile)

    if outputType == "fna":
        allSeq = gbk2fasta(inputFile)
        suffix = '.fna'
    elif outputType == "ffn":
        allSeq = gbk2ffn(inputFile)
        suffix = '.ffn'
    elif outputType == "faa":
        allSeq = gbk2faa(inputFile)
        suffix = '.faa'

    SeqIO.write(allSeq, outputFile+suffix, 'fasta')
