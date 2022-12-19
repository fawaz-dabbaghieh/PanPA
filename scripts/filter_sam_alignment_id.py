import sys

SAM_NM = 11
SAM_MD = 12
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


def get_match_miss_md(md_string):
	number = ""
	matches = []
	mismatches = 0
	for c in md_string:
		if c.isnumeric():
			number += c
		elif c == "^":  # deletion
			if number:
				matches.append(int(number))
			number = ""
			continue
		else:
			if number:
				matches.append(int(number))
			number = ""
			mismatches += 1
	if number:
		matches.append(int(number))
	return sum(matches), mismatches


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
		if not l:  # skip empty line if there's any
			continue
		if l[1] == "4":
			continue
		cigar = l[SAM_CIGAR]
		cigar_dict = cigar_to_dict(cigar)
		md = l[SAM_MD].split(":")[-1]
		matches, mismatches = get_match_miss_md(md)
		nm = int(l[SAM_NM].split(":")[-1])
		a_block_len = cigar_a_length(cigar_dict)
		#a_id = round((a_block_len - nm)/a_block_len, 5)
		a_id = (matches - nm)/a_block_len
		if a_id >= threshold:
			print(line.strip())
