import gzip
import os
import sys
import shutil


def read_fasta(fasta_file):
	sequences = dict()
	seqs = []
	seq_name = ""
	with open(os.path.join(in_dir, fasta_file), "r") as infile:
		for line in infile:
			line = line.strip()
			if not line:
				continue
	
			if line.startswith(">"):
				if len(seqs) != 0:
					sequences[seq_name] = "".join(seqs)
					seq_name = line[1:]
					seqs = []
				else:
					seq_name = line[1:]
			else:
				seqs.append(line)
		if seqs:
			sequences[seq_name] = "".join(seqs)
		return sequences
if len(sys.argv) < 3:
	print("You need to give an input and output directories")
	sys.exit()

in_dir = sys.argv[1]
out_dir = sys.argv[2]
tomove_files = []
for f in os.listdir(in_dir):
	if len(read_fasta(f)) == 1:
		tomove_files.append(os.path.join(in_dir, f))

for f in tomove_files:
	if os.path.exists(f):
		shutil.move(f, out_dir)
