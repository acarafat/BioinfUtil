'''
For a set of fasta files in a directory, retrieve number of nucleotide per fasta file.
This is helpful to get estimate of genome size.

Usage in commandline:
python get_fasta_stat.py <directory> <fasta_suffix>
'''

from Bio import SeqIO
from os import listdir
from sys import argv


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


def main():
  genome_size(argv[1], argv[2])

  
if __name__ == "__main__":
  main()
