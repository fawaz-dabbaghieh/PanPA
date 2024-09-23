import os
import sys
from PanPA.read_fasta import *


if len(sys.argv) < 2:
	print("You need to give the input directory with the fasta files")
	sys.exit()

if not sys.argv[1].endswith(os.sep):
	sys.argv[1] = sys.argv[1] + os.sep

for f in os.listdir(sys.argv[1]):
	with open("tmp", "w") as outfile:
		for seq_id, seq in read_fasta_gen(sys.argv[1] + f):
			outfile.write(f"{seq_id}\n")
			outfile.write(f"{seq.upper()}\n")

		os.rename('tmp', sys.argv[1] + f)
