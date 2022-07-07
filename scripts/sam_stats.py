import os, sys
import pickle
from collections import defaultdict
import pdb


def align_block(cigar):
	cigar = cigar.replace("M", ",").replace("I", ",").replace("D", ",").replace("S", ",").replace("H", ",")
	return sum([int(x) for x in cigar.split(",") if x])


in_sam = sys.argv[1]
"""
1139999 0  (I think aligned)
1292610 16 (aligned but reversed)
 125298 2048 (supplementary alignment)
 141454 2064 (supplementary alignment and reversed)
2407564 4 (not aligned)
"""

primary = 0
primary_rev = 16
supp = 2048
supp_rev = 2064
not_aligned = 4
stats = defaultdict(int, {0:0, 16:0, 2048:0, 2064:0, 4:0})
align_scores = set()


FLAG = 1
CIGAR = 5
NM = 11

with open(in_sam, "r") as infile:
	for l in infile:
		if l.startswith("@"):
			continue
		l = l.strip().split()

		stats[int(l[FLAG])] += 1

		if l[FLAG] == "4":
			continue

		block_size = align_block(l[CIGAR])
		nm = int(l[NM].split(":")[-1])
		a_score = round((block_size - nm)/block_size, 3)
		align_scores.add(a_score)

print(f"The input file: {in_sam}")
print(f"Number of primary alignments: {stats[primary]}")
print(f"Number of primary alignments reverse: {stats[primary_rev]}")
print(f"Number of supplementary alignments: {stats[supp]}")
print(f"Number of supplementary alignments rev: {stats[supp_rev]}")
print(f"Number of not aligned: {stats[not_aligned]}")
with open("alignment_scores_distribution_list.pickle", "bw") as outfile:
	pickle.dump(list(align_scores), outfile)
