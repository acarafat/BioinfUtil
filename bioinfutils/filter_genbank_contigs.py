from Bio import SeqIO
import argparse

def filter_genbank(input_file, output_file, contig_ids):
    """
    Filter GenBank file by contig IDs.

    Parameters:
    - input_file: Path to the input GenBank file.
    - output_file: Path to the output GenBank file.
    - contig_ids: List of contig IDs to filter by.
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # Parse the input GenBank file
        records = SeqIO.parse(infile, "genbank")

        # Filter records by contig IDs
        filtered_records = [record for record in records if record.id in contig_ids]

        # Write the filtered records to the output file
        SeqIO.write(filtered_records, outfile, "genbank")

def main():
    parser = argparse.ArgumentParser(description="Filter GenBank file based on contig IDs.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input GenBank file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output GenBank file.")
    parser.add_argument("-c", "--contigs", nargs='+', required=True, help="List of contig IDs to filter by.")

    args = parser.parse_args()

    # Call the filter function
    filter_genbank(args.input, args.output, args.contigs)

if __name__ == "__main__":
    main()
