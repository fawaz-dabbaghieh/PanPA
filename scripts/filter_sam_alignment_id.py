import os
import sys
import pickle
import pdb


SAM_NM = 11
SAM_CIGAR = 5


def cigar_to_dict(cigar):
	cigar_dict = dict()
	number = ""
	for c in cigar:
		if c.isnumeric():
			number += c
		else:
			if c not in cigar_dict:
				cigar_dict[c] = int(number)
			else:
				cigar_dict[c] += int(number)
			number = ""
	return cigar_dict


def cigar_a_length(cigar_dict, graph_aligner=False):
	length = 0
	if "I" in cigar_dict:
		length += cigar_dict["I"]
	if "D" in cigar_dict:
		length += cigar_dict["D"]
	if not graph_aligner:
		length += cigar_dict["M"]
	else:
		length += cigar_dict["="]
		try:
			length += cigar_dict["X"]
		except KeyError:
			pass  # in case there were no mismatches
	return length


if len(sys.argv) < 3:
	print("You need to give the sam file as input and the threshold 0-1")
	sys.exit()
sam_file = sys.argv[1]
threshold = float(sys.argv[2])

with open(sam_file, "r") as infile:
	for line in infile:
		if line.startswith("@"):
			print(line.strip())
			continue
		l = line.strip().split()
		cigar = l[SAM_CIGAR]
		cigar_dict = cigar_to_dict(cigar)
		nm = int(l[SAM_NM].split(":")[-1])
		a_block_len = cigar_a_length(cigar_dict)
		a_id = round((a_block_len - nm)/a_block_len, 5)
		if a_id >= threshold:
			print(line.strip())
