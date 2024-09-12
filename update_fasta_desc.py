import argparse
from Bio import SeqIO
from os import listdir


def enlist_fasta(dir_fasta, suffix):
    '''
    dir_fasta: A directory path that contains multiple fasta files. Each fasta should contain isolates for specific genes.
    suffix: Extension of the set of fasta files.
    OUTPUT: List containing the names of the fasta files in that directory.
    '''
    gene_list = []
    if not dir_fasta.endswith('/'):
        dir_fasta = dir_fasta + '/'
    for filename in listdir(dir_fasta):
        if filename.endswith(suffix):
            gene_list.append(dir_fasta + filename)
    return gene_list

def update_fasta_id(fasta_file, split_symbol, position_to_remove):
    '''
    INPUT: fasta file
    Output: fasta file containing ID without the contig ID
    '''
    seq_list = {}
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
        seq_id = seq_record.id
        # Update the fasta description
        if split_symbol in seq_id:
            seq_record.id = '_'.join(seq_id.split(split_symbol)[:-1])
            seq_record.description = ''

        # Update the dictionary
        if seq_record.id not in seq_list:
            seq_list[seq_record.id] = seq_record
        else:
            # Keep only the longest sequence if there are multiple entries for the same ID
            if len(str(seq_record.seq)) > len(str(seq_list[seq_record.id].seq)):
                seq_list[seq_record.id] = seq_record
    return seq_list


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Update FASTA descriptions based on split patterns.")

    # Positional arguments
    parser.add_argument("--input_dir", "-i", type=str, help="Path to the input directory containing FASTA files")
    parser.add_argument("--fasta_suffix", "-t",  type=str, help="Suffix of the FASTA files to be processed")

    # Optional arguments with flags
    parser.add_argument("--split_symbol", "-s", type=str, required=True, help="Symbol used to split the FASTA descriptions")
    parser.add_argument("--position_to_keep", "-p", type=int, required=True, help="Position to keep after splitting the description")
    parser.add_argument("--updated_fasta_suffix", "-u", type=str, required=True, help="Suffix for the updated FASTA output files")

    # Parse arguments
    args = parser.parse_args()

    # Get list of fasta files in the directory
    fasta_list = enlist_fasta(args.input_dir, args.fasta_suffix)

    # Process each fasta file
    for fasta_file in fasta_list:
        updated_seq_list = update_fasta_id(fasta_file, args.split_symbol, args.position_to_keep)

        # Write updated sequences to new file
        output_file = fasta_file + args.updated_fasta_suffix
        SeqIO.write(updated_seq_list.values(), output_file, 'fasta')

    print(f"FASTA files updated and saved with suffix '{args.updated_fasta_suffix}'")
