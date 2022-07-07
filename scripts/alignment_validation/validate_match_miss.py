import os
import sys
import pdb

"""
this script here should just take the gaf and looks at the g_id tag and see if the graph id
without the .gfa is present in the sequence name, if yes then it's a valid match
if not then it reports it somehow

what is important is to count the correct matches against the unique sequence names
because one sequence could align more than once, to the correct graph then the not
correct one, but if it aligns to the correct one it's valid
"""


def return_name(alignment):
	# return just the name of the sequence without the graph name attached to it
	# so can be matched later
	name = alignment.split("\t")[0]
	name = name.split(".")[0]
	name = name.split("_")
	for idx, n in enumerate(name):
		if "-" in n:
			break
	if "N" not in name[idx+1]:
		return "_".join(name[:idx + 2])
	else:
		return "_".join(name[:idx+1])


matches = set()
with open(sys.argv[1], "r") as infile:
	for l in infile:
		l = l.strip()
		matches.add(return_name(l))


mismatches = set()
with open(sys.argv[2], "r") as infile:
	for l in infile:
		l = l.strip()
		mismatches.add(return_name(l))

diff = mismatches - matches
for d in diff:
	print(d)
