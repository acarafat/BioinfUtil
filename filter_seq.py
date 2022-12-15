'''
Functions to filter files containing sequences

Update:
12/15/22 - Filtering fasta funciton updated, bug fixed. 
'''

from Bio import SeqIO

def filterFasta(seqID_list, fasta_target): 
    '''
    INPUT: 
      seqID: list of sequence ID or fasta identifier
      fasta_target: Fasta file containing sequences which need to be filtered
    OUTPUT: Fasta file containing sequences of those ID
    '''
    filtered_entry = []
    for seq_record in SeqIO.parse(fasta_target, 'fasta'):
        if seq_record.id in seqID_list:
            filtered_entry.append(seq_record)
    return filtered_entry
