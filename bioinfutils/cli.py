from program_utils import run_program

import sys

def main():
    if len(sys.argv) < 2:
        print("Usage: bioinfutil <command> [<args>]")
        sys.exit(1)

    command = sys.argv[1]

    if command == 'change_gbk_origin':
        run_program('change_gbk_origin')
        #import api.change_gbk_origin
        #api.change_gbk_origin.main()  # Call the main function 
    elif command == 'filter_fasta':
        run_program('filter_fasta')
        #import api.filter_fasta
        #api.filter_fasta.main()
    elif command == 'filter_genbank_contigs':
        run_program('filter_genbank_contigs')
        #import api.filter_genbank_contigs
        #api.filter_genbank_contigs.main()
    elif command == 'filter_var_sites':
        run_program('filter_var_sites')
        #import api.filter_var_sites
        #api.filter_var_sites.main()
    elif command == 'gbk_reverse_complement':
        run_program('gbk_reverse_complement')
        #import api.gbk_reverse_complement
        #api.gbk_reverse_complement.main()
    elif command == 'gbk2fasta':
        run_program('gbk2fasta')
        #import api.gbk2fasta
        #api.gbk2fasta.main()
    elif command == 'gene_pa_matrix':
        run_program('gene_pa_matrix')
        #import api.gene_pa_matrix
        #api.gene_pa_matrix.main()
    elif command == 'get_fasta_stat':
        run_program('get_fasta_stat')
        #import api.get_fasta_stat
        #api.get_fasta_stat.main()
    elif command == 'multi_fasta_concatenator':
        run_program('multi_fasta_concatenator')
        #import api.multi_fasta_concatenator
        #api.multi_fasta_concatenator.main()
    elif command == 'remove_from_fasta':
        run_program('remove_from_fasta')
        #import api.remove_from_fasta
        #api.remove_from_fasta.main()
    elif command == 'splti_fasta_desc':
        run_program('splti_fasta_desc')
        #import api.split_fasta_desc
        #api.split_fasta_desc.main()
    elif command == 'split_fasta_seq':
        run_program('split_fasta_seq')
        #import api.split_fasta_seq
        #api.split_fasta_seq.main()
    elif command == 'update_fasta_seqid':
        run_program('update_fasta_seqid')
        #import api.update_fasta_seqid
        #api.update_fasta_seqid.main
    else:
        print(f"Invalid command: {command}")
        sys.exit(1)

if __name__ == '__main__':
    main()