import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: bioinfutil <command> [<args>]")
        sys.exit(1)

    command = sys.argv[1]

    if command == 'change_gbk_origin':
        import change_gbk_origin
        change_gbk_origin.main()  # Call the main function 
    elif command == 'filter_fasta':
        import filter_fasta
        filter_fasta.main()
    elif command == 'filter_genbank_contigs':
        import filter_genbank_contigs
        filter_genbank_contigs.main()
    elif command == 'filter_var_sites':
        import filter_var_sites
        filter_var_sites.main()
    elif command == 'gbk_reverse_complement':
        import gbk_reverse_complement
        gbk_reverse_complement.main()
    elif command == 'gbk2fasta':
        import gbk2fasta
        gbk2fasta.main()
    elif command == 'gene_pa_matrix':
        import gene_pa_matrix
        gene_pa_matrix.main()
    elif command == 'get_fasta_stat':
        import get_fasta_stat
        get_fasta_stat.main()
    elif command == 'multi_fasta_concatenator':
        import multi_fasta_concatenator
        multi_fasta_concatenator.main()
    elif command == 'remove_from_fasta':
        import remove_from_fasta
        remove_from_fasta.main()
    elif command == 'splti_fasta_desc':
        import split_fasta_desc
        split_fasta_desc.main()
    elif command == 'split_fasta_seq':
        import split_fasta_seq
        split_fasta_seq.main()
    elif command == 'update_fasta_seqid':
        import update_fasta_seqid
        update_fasta_seqid.main
    else:
        print(f"Invalid command: {command}")
        sys.exit(1)

if __name__ == '__main__':
    main()