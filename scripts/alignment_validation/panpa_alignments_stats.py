import os
import sys
import pdb


def average(some_list):
	n = len(some_list)
	return sum(some_list)/n


in_gfa = sys.argv[1]
in_sequences= sys.argv[2]
prefix = sys.argv[3]

uniq_alignments = dict()

with open(in_gfa, "r") as infile:
	for l in infile:
		l = l.strip().split("\t")
		a_name = l[0]
		g_name = l[-1].split(":")[-1]
		g_name = "_".join(g_name.split("_")[0:4])
		g_name = g_name.split(".")[0]
		if g_name in a_name:  # match
			if g_name in uniq_alignments:
				uniq_alignments[a_name].append(True)
			else:
				uniq_alignments[a_name] = [True]

		else:  # mismatch
			if a_name not in uniq_alignments:
				uniq_alignments[a_name] = [False]
			else:
				uniq_alignments[a_name].append(False)


original_sequences = set()
with open(in_sequences, 'r') as infile:
	for l in infile:
		l = l.strip()
		if l.startswith(">"):
			original_sequences.add(l[1:])

n_matches = 0
n_mismatches = 0
n_uniq_alignments = len(uniq_alignments)
matches_dist = []
not_aligned = 0

mismatches_file = open(prefix + "_mismatched_ids.txt", "w")
not_aligned_file = open(prefix + "_not_aligned.txt", 'w')

for k, a in uniq_alignments.items():
	if k not in original_sequences:
		not_aligned += 1
		not_aligned_file.write(k + "\n")

	matches_dist.append(len(a))
	if True in a:
		n_matches += 1  # one match is all I want
	else:
		n_mismatches += 1
		mismatches_file.write(k + "\n")

mismatches_file.close()
not_aligned_file.close()

print(f"Number of matches: {n_matches}")
print(f"Number of mismatches: {n_mismatches}")
print(f"Number of unique alignments {len(uniq_alignments)}")
print(f"Average number of alignments per sequence {average(matches_dist)}")
print(f"Number of un-aligned sequences: {not_aligned}")

