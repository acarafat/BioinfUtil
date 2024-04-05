#!/usr/bin/python

import argparse, sys, os.path
from Bio import SeqIO
from Bio.SeqFeature import *

def open_sequence(arg):
	try:
		return SeqIO.read(arg, 'genbank')
	except ValueError as e:
		print('Error while parsing \'{}\': {}'.format(arg, e.message))
		sys.exit(-1)

def get_default_output(input_file, origin):
	(root, ext) = os.path.splitext(input_file)
	return "{}_{}{}".format(root, origin+1, ext)
	

def change_feature_location(f, origin, record):
	
	L = len(record.seq)

	parts = [p + (L-origin) for p in f.location.parts]

	parts = [ (FeatureLocation(p.start-L,p.end-L,strand=p.strand) 
						 if p.start > L else p) for p in parts]
	
	if len(parts) > 1:
		#see if any join
		_parts = []
		last = None
		for part in parts:
			if last:
				#if the part will join with the last one
				if last.end == part.start and last.strand == part.strand:
					last = FeatureLocation(last.start, part.end, strand=part.strand)
					continue
				else:
					_parts.append(last)
			last = part
		_parts.append(last)


		parts = _parts

	#split parts which are past the end
	parts = [(CompoundLocation([FeatureLocation(p.start,L,strand=p.strand),
														 FeatureLocation(0,p.end-L,strand=p.strand)]) 
						if p.end > L else p) for p in parts]

	if len(parts) == 1:
		f.location = parts[0]
	else:
		f.location = CompoundLocation(parts)

	return f

def main():
	parser = argparse.ArgumentParser(description='Change the origin of a circular	DNA genbank file')
	parser.add_argument('origin', 
											metavar='N', 
											type=int, 
											help="The 1-offset location of the new origin")
	parser.add_argument('input', 
											type = str,
											help = 'The sequence file')
	parser.add_argument('output',
											help = 'File to write to',
											nargs = '?',
											default = None)

	args = parser.parse_args()

	record = open_sequence(args.input)
	origin = args.origin -1
	output = args.output if args.output else get_default_output(args.input,	origin)

	if origin > len(record.seq):
		print('Error: New origin is larger than sequence length')
		sys.exit(-1)

	
	#move sequence
	record.seq = record.seq[origin:] + record.seq[:origin]

	#move features
	record.features = [change_feature_location(f, origin, record) 
										 for f in record.features]

	SeqIO.write(record, output, 'genbank')


if __name__=='__main__':
	main()
