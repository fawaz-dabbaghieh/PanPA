import sys
import os
from PanPA.read_fasta import *
import pdb

if len(sys.argv) < 3:
	print("You need to give the input reference fasta file and length of seq to extract")
	sys.exit()

seqs = read_fasta(sys.argv[1])

ref = seqs[list(seqs.keys())[0]]
assert len(ref) > 1000_000
start = 0


while True:
	# pdb.set_trace()
	try:
		assert start + int((int(sys.argv[2]))/2) < len(ref)
		to_print = ref[start:start + int(sys.argv[2])]
	except AssertionError:
		print(f">ref_{start}_{len(ref)}")
		print(ref[start:])
		break
	print(f">ref_{start}_{start + int(sys.argv[2])}")
	print(to_print)
	start += int((int(sys.argv[2]))/2)