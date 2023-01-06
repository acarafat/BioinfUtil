'''
Functions to filter files containing sequences

Update:
01/05/23 - New function added to only check duplicated entry and filter
12/16/22 - Filtering fasta funciton updated to check for duplicated entry
12/15/22 - Filtering fasta funciton updated, bug fixed. 
'''

from Bio import SeqIO

def byList(seqID_list, fasta_target):
    
    '''
    INPUT: 
      seqID: list of sequence ID or fasta identifier
      fasta_target: Fasta file containing sequences which need to be filtered
    OUTPUT: Fasta file containing sequences of those ID
    '''
    filtered_entry = {}
    for seq_record in SeqIO.parse(fasta_target, 'fasta'):
        if seq_record.id in seqID_list:
            if seq_record.id not in filtered_entry.keys():
                filtered_entry[seq_record.id] = seq_record
            else:
                print('duplication: ', seq_record.id)
    return list(filtered_entry.values())

def bySize(fasta_target, new_fasta):
    '''
    INPUT: Fasta file containing duplicate sequences of multiple lenght
    ALOGRITHM: Only keep the longest sequences
    Output: Fasta file contiaing unique sequences
    '''
    unique_seqs = {}
    for seq_record in SeqIO.parse(fasta_target, 'fasta'):
        if seq_record.id not in unique_seqs.keys():
            unique_seqs[seq_record.id] = seq_record
        else:
            if len(str(seq_record.seq)) > len(str(unique_seqs[seq_record.id])):
                unique_seqs[seq_record.id] = seq_record
    SeqIO.write(unique_seqs.values(), new_fasta, 'fasta')
    pass
            
