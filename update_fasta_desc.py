from Bio import SeqIO
from sys import argv

'''
What it does:
- Remove contig number from seqID
- Only keep longest sequence if there are multiple entry for same seqID

How to use this script in command line:
python update_fasta_desc.py <fasta_file_input> <fasta_file_output>
'''


def update_fasta_id(fasta_file):
        '''
        INPUT: fasta file
        Output: fasta file contianing id without the contig id
        '''
        seq_list = {}
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
                seq_id = seq_record.id
                if '_' in seq_id:
                        seq_record.id = seq_id.split('_', -1)[0]

                if seq_record.id not in seq_list.keys():
                        seq_list[seq_record.id] = seq_record
                else:
                        if len(str(seq_record.seq)) > len(str(seq_list[seq_record.id].seq)):
                                seq_list[seq_record.id] = seq_record
        return seq_list


if __name__ == "__main__":
        updated_seq_list = update_fasta_id(argv[1])
        SeqIO.write(updated_seq_list.values(), argv[2], 'fasta')
