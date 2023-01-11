'''
What it does:
- Use provided pattern to split fasta description
- Remove contig number from seqID
- Only keep longest sequence if there are multiple entry for same seqID


How to use this script in command line:
python update_fasta_desc.py <fasta_file_input> <split_symbol> <position_to_keep> <fasta_file_output>

For example, let's say a fasta description is the following:

>05LoS16R10_36_00608__05LoS16R10_36

Here the split symbol will be '__' and position to keep will be 1 (count started from 0)
'''

from Bio import SeqIO
from sys import argv


def update_fasta_id(fasta_file, split_symbol, position_to_keep):
        '''
        INPUT: fasta file
        Output: fasta file contianing id without the contig id
        '''
        seq_list = {}
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
                seq_id = seq_record.id
                if split_symbol in seq_id:
                        seq_record.id = seq_id.split(split_symbol)[position_to_keep]

                if seq_record.id not in seq_list.keys():
                        seq_list[seq_record.id] = seq_record
                else:
                        if len(str(seq_record.seq)) > len(str(seq_list[seq_record.id].seq)):
                                seq_list[seq_record.id] = seq_record
        return seq_list


if __name__ == "__main__":
        updated_seq_list = update_fasta_id(argv[1], argv[2], argv[3])
        SeqIO.write(updated_seq_list.values(), argv[4], 'fasta')
