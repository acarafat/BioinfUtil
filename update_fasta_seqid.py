from Bio import SeqIO
import argparse


def update_fasta_ids(input_file, output_file):
  """
  Updates FASTA sequence IDs in a file.

  Args:
    input_file: Path to the input FASTA file.
    output_file: Path to the output FASTA file.
  """
  # Get filename without extension
  filename = input_file.split(".")[0]

  with open(input_file) as handle, open(output_file, "w") as output:
    all_seq_record = []
    for record in SeqIO.parse(handle, "fasta"):
      new_id = f"{filename}__{record.id}"
      record.id = new_id
      all_seq_record.append(record)
    SeqIO.write(all_seq_record, output, "fasta")


if __name__ == "__main__":
  # Define arguments
  parser = argparse.ArgumentParser(description="Update FASTA sequence IDs")
  parser.add_argument( "-i", "--input", help="Path to the input FASTA file")
  parser.add_argument("-o",  "--output", help="Path to the output FASTA file")
  args = parser.parse_args()

  # Update FASTA file
  update_fasta_ids(args.input, args.output)

  print(f"FASTA sequence IDs updated in: {args.output}")
