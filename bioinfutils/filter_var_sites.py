from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import argparse

import os
import glob

def detect_variable_sites(aa_aln):
    '''
    INPUT: Amino acid alignment file in fasta format
    OUTPUT: Position of non-synonymous sites
    '''
    # Read the alignment file
    alignment = AlignIO.read(aa_aln, 'fasta')

    # Number of sequences and their length
    num_sequences = len(alignment)
    seq_length = alignment.get_alignment_length()

    variable_sites = []

    # Iterate through each column (site) in the alignment
    for i in range(seq_length):
        # Get the set of amino acids at this position
        amino_acids = set(record.seq[i] for record in alignment)

        # If there's more than one unique amino acid, it's a variable site
        if len(amino_acids) > 1:
            variable_sites.append(i) # positiona are 0-based

    return variable_sites

def aa_to_na_pos(variable_sites):
    '''
    INPUT: List of AA position(s)
    OUTPUT: List of NA codon positions
    '''
    codon_positions = []
    for i in variable_sites:
        n = i*3
        codon_positions.extend([n, n+1, n+2])

    return codon_positions


def remove_variable_sites(alignment, pos_list):
    '''
    INPUT: List of position and alignment file in fasta format
    OUTPUT: A new alignment file where positions from pos_list were removed
    '''

    # Read alignment file
    alignment = AlignIO.read(alignment, "fasta")

    # Create a new modified alignment object
    modified_alignment = MultipleSeqAlignment(
            [
                SeqRecord(
                    Seq("".join(base for i, base in enumerate(record.seq) if i not in pos_list)),
                    id = record.id,
                    description=record.description
                )
                for record in alignment
            ]
        )

    return modified_alignment


def remove_nonsynonymous_sites(codon_alignment, translated_alignment, output):
    '''
    INPUT: A codon alignment file and a tranlated alignment file in fasta format.
    OUTPUT: A nucleotide alignment file where all non-synonymous sites has been removed
    '''
    variable_sites = detect_variable_sites(translated_alignment)
    codon_sites = aa_to_na_pos(variable_sites)
    modified_aln = remove_variable_sites(translated_alignment, codon_sites)
    with open(output, 'w') as output_handle:
        AlignIO.write(modified_aln, output_handle, 'fasta')

    pass


def main(args=None):
    parser = argparse.ArgumentParser(description="Detect variable sites in an amino acide alignment.")

    parser.add_argument("--input_aa_fasta", type=str, help="Path to amino acid alignment file in fasta format")
    parser.add_argument("--input_na_fasta", type=str, help="Path to nucleotide acid codon alignment file in fasta format")
    parser.add_argument("--output_fasta", type=str, help="Path to modified codon alignment output file")

    parser.add_argument("--input_aa_dir", type=str, help="Input directory path containing translated codon alignment fasta")
    parser.add_argument("--input_na_dir", type=str, help="Input directory path containing codon alignment fasta file")
    parser.add_argument("--input_aa_suffix", type=str, help="Common suffix of AA fasta file input")
    parser.add_argument("--input_na_suffix", type=str, help="Common suffix of NA fasta file input")

    args = parser.parse_args(args)

    # Case: only AA and NA files are provided for a single strain
    if args.input_aa_fasta != None:
        remove_nonsynonymous_sites(args.input_na_fasta, args.input_aa_fasta, args.output_fasa)


    # Case: Directories to AA and NA files are provided for multiple strains/genes
    else:
        aa_files = glob.glob(os.path.join(args.input_aa_dir, "*"+args.input_aa_suffix))
        na_files = glob.glob(os.path.join(args.input_na_dir, "*"+args.input_na_suffix))

        # Pair the codon and translated alignment file, matching by file name prefix)
        fasta_pairs = {}

        for i in aa_files:
            basename = i.split(args.input_aa_suffix)[0]
            if '/' in basename:
                basename = basename.split('/')[-1]
            for j in na_files:
                if basename in j:
                    fasta_pairs[basename] = [i, j]
                    break

        for k in fasta_pairs.keys():
            remove_nonsynonymous_sites(fasta_pairs[k][1], fasta_pairs[k][0], k+'.neutral.fna')


if __name__ == "__main__":
    main()