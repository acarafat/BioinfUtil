import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: bioinfutil <command> [<args>]")
        sys.exit(1)

    command = sys.argv[1]

    if command == 'change_gbk_origin':
        from . import change_gbk_origin
        change_gbk_origin.main()  # Call the main function 
    elif command == 'filter_fasta':
        from . import filter_fasta
        filter_fasta.main()
    elif command == 'filter_genbank_contigs':
        from . import filter_genbank_contigs
        filter_genbank_contigs.main()
    elif command == 'filter_var_sites':
        from . import filter_var_sites
        filter_var_sites.main()
    elif command == 'gbk_reverse_complement':
        from . import gbk_reverse_complement
        gbk_reverse_complement.main()
    elif command == 'gbk2fasta':
        from . import gbk2fasta
        gbk2fasta.main()
    elif command == 'gene_pa_matrix':
        from . import gene_pa_matrix
        gene_pa_matrix.main()
    elif command == 'get_fasta_stat':
        from . import get_fasta_stat
        get_fasta_stat.main()
    elif command == 'multi_fasta_concatenator':
        from . import multi_fasta_concatenator
        multi_fasta_concatenator.main()
    elif command == 'remove_from_fasta':
        from . import remove_from_fasta
        remove_from_fasta.main()
    elif command == 'split_fasta_desc':
        from . import split_fasta_desc
        split_fasta_desc.main()
    elif command == 'split_fasta_seq':
        from . import split_fasta_seq
        split_fasta_seq.main()
    elif command == 'update_fasta_seqid':
        from . import update_fasta_seqid
        update_fasta_seqid.main()
    else:
        print(f"Invalid command: {command}")
        sys.exit(1)

if __name__ == '__main__':
    main()