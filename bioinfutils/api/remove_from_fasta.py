from Bio import SeqIO
import argparse

def remove_from_fasta(fasta_file, sequences_to_remove, output_file):
  """
  This function removes sequences from a FASTA file based on a list of sequence IDs.

  Args:
    fasta_file: Path to the input FASTA file.
    sequences_to_remove: List of sequence IDs to remove.
    output_file: Path to the output FASTA file.
  """
  # Read the FASTA file
  with open(fasta_file, "r") as handle:
    records = list(SeqIO.parse(handle, "fasta"))

  # Read sequence IDs to remove from a file
  with open(sequences_to_remove, "r") as handle:
    sequences_to_remove = [line.strip() for line in handle.readlines()]

  # Filter records to keep only those not in the removal list
  filtered_records = [record for record in records if record.id not in sequences_to_remove]

  # Write filtered records to a new FASTA file
  with open(output_file, "w") as handle:
    SeqIO.write(filtered_records, handle, "fasta")

def main():
  # Create argument parser
  parser = argparse.ArgumentParser(description="Remove sequences from a FASTA file.")

  # Required arguments
  parser.add_argument("fasta_file", help="Path to the input FASTA file.")
  parser.add_argument("sequences_to_remove", help="Path to a file containing list of sequence IDs to remove (one ID per line).")
  parser.add_argument("output_file", help="Path to the output FASTA file.")

  # Parse arguments
  args = parser.parse_args()

  # Call function with parsed arguments
  remove_from_fasta(args.fasta_file, args.sequences_to_remove, args.output_file)

  print(f"Sequences removed from {args.fasta_file} and written to {args.output_file}")

if __name__ == "__main__":
  main()
