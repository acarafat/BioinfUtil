# Split a set of sequences in a fasta file into individual fasta files

from Bio import SeqIO
import sys

def splitFastaSeq(inputFile, outputDir, prefix=''):
    """
    Input: A fasta file containing set of sequences
    Output: Set of fasta files containing individual sequences in the outputDir with `prefix`
    """
    if not outputDir.endswith('/'):
        outputDir = outputDir + '/'
    for seq_record in SeqIO.parse(inputFile, 'fasta'):
        SeqIO.write(seq_record, outputDir+prefix+str(seq_record.id), 'fasta')

    pass

if __name__ == "__main__":
    inputFile = sys.argv[1]
    outputDir = sys.argv[2]
    prefix = sys.argv[3]
    splitFastaSeq(inputFile, outputDir, prefix)
    print('Done.')
