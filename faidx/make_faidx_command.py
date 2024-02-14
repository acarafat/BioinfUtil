help = '''
Given a coordinate file, this will genarate samtools faidx command to extract target sequences.

python3 make_faidx_command.py coordinate.txt genomes.fasta


coordinate.txt is a tab seperated file and look like this:
strain1  58700  59980
strain2  115021  113741
strain3  160429  161709
<file>  <start position>  <end position>
...  ...  ...
'''

from sys import argv


def main(position_file, target_fasta):
   '''
   INPUT: position file for samtools faidx
   OUTPUT: command file for running samtools faidx
   '''
   unique = {}
   command = ''

   with open(position_file) as handle:
      raw = handle.readlines()

   for i in raw:
      seq_id, m, n = i.split()
      genome_id = seq_id.split('_')[0]
      if genome_id not in unique.keys():
            unique[genome_id] = [seq_id, m, n, abs(int(m)-int(n))]
      else:
         if abs(int(m)-int(n)) > unique[genome_id][3]:
            unique[genome_id] = [seq_id, m, n, abs(int(m)-int(n))]
      if int(unique[genome_id][1]) > int(unique[genome_id][2]):
         unique[genome_id] = [seq_id, n, m, abs(int(m)-int(n))]

   for k in unique.keys():
      command += f"samtools faidx {target_fasta} {unique[k][0]}:{unique[k][1]}-{unique[k][2]}\n"

   with open('faidx_command.sh','w') as handle:
      handle.write(command)

   pass


if __name__ == "__main__":
  try:
    main(argv[1], argv[2])
  except:
    print(help)
