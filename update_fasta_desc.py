'''
What it does:
- Use provided pattern to split fasta description
- Keep the expected substring from provided position
- Only keep longest sequence if there are multiple entry for same seqID


How to use this script in command line:
python update_fasta_desc.py <~/path/to/input/fasta> <fasta_suffix> <split_symbol> <position_to_keep> <updated_fasta_suffix>

For example, let's say a fasta description is the following:

>05LoS16R10_36_00608__05LoS16R10_36

Here the split symbol will be '__' and position to keep will be 1 (count started from 0)
'''

from Bio import SeqIO
from sys import argv
from os import listdir


def enlist_fasta(dir_fasta, suffix):
    '''
    dir_fasta: A directory path that contains multiple fasta files. Each fasta should contains isolates for a specific genes.
    SUFFIX: Extension of the set of fasta file.
    OUTPUT: List containing name of the fasta files in that directory.
    '''
    gene_list = []
    if not dir_fasta.endswith('/'):
        dir_fasta = dir_fasta + '/'
    for filename in listdir(dir_fasta):
        if filename.endswith(suffix):
            gene_list.append(dir_fasta+filename)
    return gene_list

def update_fasta_id(fasta_file, split_symbol, position_to_keep):
        '''
        INPUT: fasta file
        Output: fasta file contianing id without the contig id
        '''
        seq_list = {}
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
                seq_id = seq_record.id
                # Update the fasta description
                if split_symbol in seq_id:
                        seq_record.id = seq_id.split(split_symbol)[int(position_to_keep)]
                        seq_record.description = ''
                        
                # Update the dictionary        
                if seq_record.id not in seq_list.keys():
                        seq_list[seq_record.id] = seq_record
                else:
                        if len(str(seq_record.seq)) > len(str(seq_list[seq_record.id].seq)):
                                seq_list[seq_record.id] = seq_record
        return seq_list


if __name__ == "__main__":
        fasta_list = enlist_fasta(argv[1], argv[2])
        for fasta_file in fasta_list:
                updated_seq_list = update_fasta_id(fasta_file, argv[3], argv[4])
                SeqIO.write(updated_seq_list.values(), fasta_file+argv[5], 'fasta')
