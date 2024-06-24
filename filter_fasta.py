'''
'''

from Bio import SeqIO
import argparse


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
            

def filter_fasta_by_id(fasta_file, id_file, output_file):
  """Filters a FASTA file by sequence IDs in another file using Biopython.

  Args:
    fasta_file: Path to the FASTA file.
    id_file: Path to the text file containing IDs (one per line).
    output_file: Path to the output FASTA file.
  """

  with open(id_file, 'r') as ids:
    desired_ids = set(line.strip() for line in ids)

  # Use Biopython's SeqIO to iterate through FASTA records
  with open(fasta_file, 'r') as fasta, open(output_file, 'w') as out:
    for record in SeqIO.parse(fasta, 'fasta'):
      if record.id in desired_ids:
        SeqIO.write(record, out, 'fasta')


if __name__ == '__main__':
  # Define argument parser (same as before)
  parser = argparse.ArgumentParser(description='Filter FASTA file by sequence IDs, size, or list')
  parser.add_argument('-f', '--fasta', required = True, help='Path to the FASTA file')
  parser.add_argument('-i', '--ids', required = False, help='Path to the file containing IDs')
  parser.add_argument('-o', '--output', required = True, help='Path to the output FASTA file')
  parser.add_argument('-n', '--option', required = True, help='Type of filter: 1 for by id in list, 2 for by size, 3 $
  parser.add_argument('-l', '--list', required =  False, help='List of contig id\'s to be filtered')

  # Parse arguments (same as before)
  args = parser.parse_args()

  # Call the filtering function
  if args.option == '3':
    filter_fasta_by_id(args.fasta, args.ids, args.output)
  elif args.option == '2':
    bySize(args.fasta, args.output)
  else:
    byList(args.list, args.fasta, args.output)

  print(f"Filtered FASTA written to: {args.output}")

