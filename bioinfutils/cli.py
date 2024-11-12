import sys
import importlib

# Mapping of subcommands to their respective module paths
commands = {
    'change_gbk_origin': 'bioinfutils.change_gbk_origin',
    'filter_fasta': 'bioinfutils.filter_fasta',
    'filter_genbank_contigs': 'bioinfutils.filter_genbank_contigs',
    'filter_var_sites': 'bioinfutils.filter_var_sites',
    'find_gbk': 'bioinfutils.find_gbk',
    'gbk_reverse_complement': 'gbk_reverse_complement',
    'gbk2fasta': 'bioinfutils.gbk2fasta',
    'gene_pa_matrix': 'bioinfutils.gene_pa_matrix',
    'get_fasta_stat': 'bioinfutils.get_fasta_stat',
    'multi_fasta_concatenator': 'bioinfutils.multi_fasta_concatenator',
    'remove_from_fasta': 'bioinfutils.remove_from_fasta',
    'split_fasta_desc': 'bioinfutils.split_fasta_desc',
    'split_fasta_seq': 'bioinfutils.split_fasta_seq',
    'update_fasta_seqid': 'bioinfutils.update_fasta_seqid'
    # Add more commands here
}

def main():
    if len(sys.argv) < 2:
        print("Usage: bioinfutils <command> [<args>]")
        sys.exit(1)

    command = sys.argv[1]
    if command in commands:
        # Import the appropriate module based on the command
        module = importlib.import_module(commands[command])
        
        # Call the main() function of the module, passing all remaining sys.argv arguments
        module.main(sys.argv[2:])  # Pass remaining arguments to module's main()
    else:
        print(f"Error: Unrecognized command '{command}'")
        print("Available commands:", ", ".join(commands.keys()))
        sys.exit(1)

if __name__ == "__main__":
    main()