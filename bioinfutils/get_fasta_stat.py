from Bio import SeqIO
from os import listdir
import argparse


def genome_size(fasta_dir, fasta_suffix):
  '''
  INPUT: A fasta file containing directory and suffix of fasta file (i.e. .fas)
  OUTPUT: A file containing genome size in each fasta file
  '''
  fasta_size = {}
  
  for fasta_file in listdir(fasta_dir):
    if fasta_file.endswith(fasta_suffix):
      for seq_record in SeqIO.parse(fasta_file, 'fasta'):
        if fasta_file not in fasta_size.keys():
         fasta_size[fasta_file] = int(len(str(seq_record.seq)))
        else:
          fasta_size[fasta_file] += int(len(str(seq_record.seq)))
    
  with open('fasta_stat.csv', 'w') as f:
    for key in fasta_size.keys():
      f.write(f"{key}, {fasta_size[key]}\n")
  pass


def main(args=None):
  parser = argparse.ArgumentParser(description='Get fasta stat')
  parser.add_argument('--input', '-i',  help='Input fasta directory.')
  parser.add_argument('--suffix', '-s', help='fasta suffix')
  args = parser.parse_args(args)

  genome_size(args.input, args.suffix)

  
if __name__ == "__main__":
  main()
