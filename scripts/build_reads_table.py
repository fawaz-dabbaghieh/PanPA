import os
import sys
import pickle
import pdb


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


def bwa_a_len(cigar_dict, line):
	"""
	To get the proper alignment length percentage, I need to see how long the sequence is
	and consider hard and soft clipping, then I can measure how much was actually aligned
	"""
	original_seq_len = 0
	aligned = 0
	seq_in_sam = len(line[9])
	if "H" in cigar_dict:  # there's a hard clip, so the original seq len is seq_in_sam + cigar_dict["H"]
		original_seq_len = seq_in_sam + cigar_dict["H"]
		if "S" in cigar_dict:  # there's soft clipping that need to be removed
			aligned = original_seq_len - cigar_dict["H"] - cigar_dict["S"]
		else:
			aligned = original_seq_len - cigar_dict["H"]
	else:
		original_seq_len = seq_in_sam  # nothing clipped from the sequence in sam
		if "S" in cigar_dict:
			aligned = original_seq_len - cigar_dict["S"]
		else:  # it is all aligned
			aligned = original_seq_len
	return round(aligned/original_seq_len, 5)
#################### table columns idx constants
A_BWA = 0
NUMB_BWA = 1
LEN_BWA = 2
ID_BWA = 3

A_GA = 4
NUMB_GA = 5
LEN_GA = 6
ID_GA = 7

A_PANPA = 8
NUMB_PANPA = 9
LEN_PANPA = 10
ID_PANPA = 11


#################### files constants
SAM_CIGAR = 5
SAM_NM = 11

GA_CIGAR = -1
GA_ID = 15

PANPA_CIGAR = 15
PANPA_ID = 14

SEQ_LEN = 1
A_START = 2
A_END = 3
####################

if len(sys.argv) < 5:
	print("You need to give <original_sequences_aligned.txt> <bwa.sam> <ga.gaf> <panpa.gaf>")
	sys.exit()
sequences_file = sys.argv[1]
bwa_file = sys.argv[2]
ga_file = sys.argv[3]
panpa_file = sys.argv[4]

big_table = dict()
big_table['header'] = ["in_bwa", "n_alignments_bwa", "bwa_alignment_length", "bwa_alignment_id",
"in_ga", "n_alignments_ga", "ga_alignment_length", "ga_alignment_id",
"in_panpa", "n_alignments_panpa", "panpa_alignment_length", "panpa_alignment_id"]

# get all seq names
with open(sequences_file, "r") as infile:
	for l in infile:
		big_table[l.strip()] = [0] * 12

# getting sam info
# Need to add a column for SAM about which flag the alignment got, so if it's aligned or not, so I can do a Venn diagram or upset plot, then plot the scores for only the aligned ones
with open(bwa_file, "r") as infile:
	for l in infile:
		if l.startswith("@"):
			continue
		l = l.strip().split()
		name = l[0]
		cigar = l[SAM_CIGAR]
		cigar_dict = cigar_to_dict(cigar)
		nm = int(l[SAM_NM].split(":")[-1])
		a_len = bwa_a_len(cigar_dict, l)
		big_table[name][A_BWA] = 1  # this sequence has been aligned by BWA
		big_table[name][NUMB_BWA] += 1  # how many time this sqeuence was aligned
		a_block_len = cigar_a_length(cigar_dict)
		a_id = round((a_block_len - nm)/a_block_len, 5)
		if big_table[name][ID_BWA] < a_id:
			big_table[name][ID_BWA] = a_id
			big_table[name][LEN_BWA] = a_len
			# big_table[name][SAM_NM] = nm
			# big_table[name][SAM_CIGAR] = cigar
		# if big_table[name][LEN_BWA] < a_len:  # only keeping the score of the better alignemnt
		# 	big_table[name][LEN_BWA] = a_len
		# 	# I think this is wrong, a_id should be (a_len - nm)/a_len

		# if big_table[name][ID_BWA] < a_id:  # only keeping the better id score
		# 	big_table[name][ID_BWA] = a_id


# getting ga gaf info
with open(ga_file, "r") as infile:
	for l in infile:
		l = l.strip().split()
		name = l[0]
		big_table[name][A_GA] = 1  # alignment exists
		# cigar = l[GA_CIGAR].split(":")[-1]
		# cigar_dict = cigar_to_dict(cigar)
		# the newer version of graph aligner have = instead of M for matches and X for mismatches
		a_len = round((int(l[A_END]) - int(l[A_START]))/int(l[SEQ_LEN]),5)
		big_table[name][NUMB_GA] += 1  # if more than once aligned
		a_id = round(float(l[GA_ID].split(":")[-1]), 5)
		if big_table[name][ID_GA] < a_id:
			big_table[name][ID_GA] = a_id
			big_table[name][LEN_GA] = a_len

		# if big_table[name][LEN_GA] < a_len:
		# 	big_table[name][LEN_GA] = a_len

		# if big_table[name][ID_GA] < a_id:
		# 	big_table[name][ID_GA] = a_id


# getting ga gaf info
with open(panpa_file, "r") as infile:
	for l in infile:
		l = l.strip().split()
		name = l[0]
		# cigar = l[PANPA_CIGAR].split(":")[-1]
		# cigar_dict = cigar_to_dict(cigar)
		# a_len = cigar_a_length(cigar_dict)
		a_len = round((int(l[A_END]) - int(l[A_START]))/int(l[SEQ_LEN]),5)
		big_table[name][A_PANPA] = 1
		big_table[name][NUMB_PANPA] += 1
		a_id = round(float(l[PANPA_ID].split(":")[-1]), 5)
		if big_table[name][ID_PANPA] < a_id:
			big_table[name][ID_PANPA] = a_id
			big_table[name][LEN_PANPA] = a_len
		# if big_table[name][LEN_PANPA] < a_len:
		# 	big_table[name][LEN_PANPA] = a_len

		# if big_table[name][ID_PANPA] < a_id:
		# 	big_table[name][ID_PANPA] = a_id

with open("big_table.pickle", "wb") as outfile:
	pickle.dump(big_table, outfile)