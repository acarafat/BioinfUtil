from Bio import SeqIO
import argparse


def split_fasta_seq(input_file, output_dir, prefix=""):
  """
  Splits a set of sequences in a FASTA file into individual FASTA files.

  Args:
    input_file: Path to the input FASTA file.
    output_dir: Path to the output directory.
    prefix: Optional prefix to add to the output filenames (default: "").
  """
  if not output_dir.endswith("/"):
    output_dir += "/"

  for seq_record in SeqIO.parse(input_file, "fasta"):
    output_file = f"{output_dir}{prefix}{seq_record.id}.fasta"
    SeqIO.write(seq_record, output_file, "fasta")


if __name__ == "__main__":
  # Define arguments
  parser = argparse.ArgumentParser(description="Split FASTA file into individual files")
  parser.add_argument("-i", "--input",  help="Path to the input FASTA file", required=True)
  parser.add_argument("-o", "--output", help="Path to the output directory", required=True)
  parser.add_argument("-p", "--prefix", help="Optional prefix for output filenames", required=False, default="")

  args = parser.parse_args()

  # Split FASTA file
  split_fasta_seq(args.input, args.output, args.prefix)

  print("Done. Individual FASTA files created in the output directory.")
