#!/bin/python

# ./inch_assem_filt.py -f [path_to_file] -s [length filter]
# This filters inchworm assembly contigs (or any fasta in 2 line format) less than X basepairs (limit is modifiable)
# Program takes input, 's', which defines the filter length below which sequences will be discarded, and 'p' the path to the input file

import argparse as ap

def get_arguments():
	parser = ap.ArgumentParser(description="Fasta/Inchworm contig size filterer")
	parser.add_argument("-p", help="directory containing the input file (use absolute path)", required=True, type=str)
	parser.add_argument("-s", help="Sequence length filter size (below which seqs are discarded)", required=True, type=int)
	return parser.parse_args()

args = get_arguments()	

def filename_creation(directory):
	filename = directory.strip('\n') + '/inchworm.DS.fa'
	outfile = directory.strip('\n') + '/inchworm.K25.L25.DS.filtered.fa'
	return filename, outfile

def contig_filter(filename, outfile, filter):

	with open(filename, 'r') as in_assem:
		out = open(outfile, 'a')
		line = in_assem.readline()
		while line:
			length = int((line.split(':')[4].strip(' ')).strip('\n'))
			if length > filter:
				out.write(line)
				line = in_assem.readline()
				while line.startswith('>') != True:
					out.write(line)
					line = in_assem.readline()

			else:
				line = in_assem.readline()
				while line.startswith('>') != True:
					if line == '':
						break
					else:
						line = in_assem.readline()

		out.close()


directory = args.p
flen = args.s
file = filename_creation(directory)[0]
ofile = filename_creation(directory)[1]
print("filtering assembly ...", file)
contig_filter(file, ofile, flen)
print("Filtering Complete")
