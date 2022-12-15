'''
Functions to filter files containing sequences
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
    if seq_record.id in seqIO_list:
      filtered_entry.append(seq_record)
      
   return filtered_entry
  
