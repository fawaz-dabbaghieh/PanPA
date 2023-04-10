import os
import sys
import pdb

translation_table = {
	'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
	'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
	'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
	'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
	'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
	'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
	'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
	'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
	'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
	'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
	'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
	'TAC': 'Y', 'TAT': 'Y', 'TAA': '[', 'TAG': '[',
	'TGC': 'C', 'TGT': 'C', 'TGA': '[', 'TGG': 'W',
}

def print_dp(dp_table):
	for i in range(len(dp_table)):
		print(" ".join([str(x) for x in dp_table[i]]))


def fs_aware_sw_codon_start(seq, target, match=2, mismatch=-2, gap=-1, fs=-2):
	dp_table = []

	max_score = 0
	max_score_loc = []

	for i in range(len(seq) + 4):
		dp_table.append([0] * (len(target) + 1))

	for i in range(len(seq) - 2):  # rows
		
		seq_char = translation_table["".join([seq[i], seq[i+1], seq[i+2]])]
		i = i + 4  # skipping the initialization lines

		for j in range(len(target)):  # columns
			target_char = target[j]
			j = j + 1  # skipping the initialization column

			if seq_char == target_char:
				score = match
			else:
				score = mismatch

			dp_table[i][j] = max(
				# same frame
				dp_table[i-3][j-1] + score,
				dp_table[i][j-1] + gap,
				dp_table[i-3][j] + gap,

				# other frames
				dp_table[i-4][j-1] + fs,
				dp_table[i-2][j-1] + fs,
				0
				)

			# updating the max score and max score location
			if dp_table[i][j] > max_score:
				max_score = dp_table[i][j]
				max_score_loc = [(i, j)]
			elif dp_table[i][j] == max_score:
				max_score_loc.append((i,j))


	# print(max_score, max_score_loc)

	# traceback from the locations of max score and find out where I went
	for start_i, start_j in max_score_loc:
		out_seq = ""
		out_target = ""
		diag_type = 3  # to tell me how many letters to add to out_seq when jumping

		while start_i > 3 or start_j > 0 or dp_table[start_i][start_j] != 0:
			# pdb.set_trace()
			i = start_i - 4  # because of the 4 initialization lines
			j = start_j - 1  # because of the 1 initialization column

			seq_char = translation_table["".join([seq[i], seq[i+1], seq[i+2]])]
			target_char = target[j]

			if seq_char == target_char:
				score = match
			else:
				score = mismatch

			######## within frame operations
			if dp_table[start_i][start_j] == dp_table[start_i-3][start_j-1] + score:  # match
				# I need to know what was my previous jump to know which letters to add to output
				if diag_type == 3:
					out_seq += seq[i+2] + seq[i+1] + seq[i]
				elif diag_type == 4:
					out_seq += seq[i+3] + seq[i+2] + seq[i+1] + seq[i]
				else:
					out_seq += seq[i+1] + seq[i]

				out_target += target_char
				diag_type = 3
				start_i, start_j = start_i-3, start_j-1


			elif dp_table[start_i][start_j] == dp_table[start_i][start_j-1] + score:  # deletion in seq
				out_seq += "---"
				out_target += target_char
				start_i, start_j = start_i, start_j-1


			elif dp_table[start_i][start_j] == dp_table[start_i-3][start_j] + score:  # insertion in seq
				out_seq += seq[i+2] + seq[i+1] + seq[i]
				out_target += "-"
				start_i, start_j = start_i-3, start_j


			####### frameshift operations
			elif dp_table[start_i][start_j] == dp_table[start_i-4][start_j-1] + fs:  # insertion shift, -1 frame
				if diag_type == 3:
					out_seq += seq[i+2] + seq[i+1] + seq[i]
				elif diag_type == 4:
					out_seq += seq[i+3] + seq[i+2] + seq[i+1] + seq[i]
				else:
					out_seq += seq[i+1] + seq[i]

				diag_type = 4
				out_target += target_char
				start_i, start_j = start_i-4, start_j-1


			elif dp_table[start_i][start_j] == dp_table[start_i-2][start_j-1] + fs:  # deletion shift, +1 frame
				if diag_type == 3:
					out_seq += seq[i+2] + seq[i+1] + seq[i]
				elif diag_type == 4:
					out_seq += seq[i+3] + seq[i+2] + seq[i+1] + seq[i]
				else:
					out_seq += seq[i+1] + seq[i]
				diag_type = 2
				out_target += target_char

				start_i, start_j = start_i-2, start_j-1

		print(out_seq[::-1])
		print(out_target[::-1])

	return dp_table

def fs_aware_sw_codon_end(seq, target, match=2, mismatch=-2, gap=-1, fs=-2):
	dp_table = []

	max_score = 0
	max_score_loc = []

	for i in range(len(seq) + 2):
		dp_table.append([0] * (len(target) + 1))

	for i in range(len(seq) - 2):  # skipping two because 2 characters are not part of a codon

		i = i + 2  # skipping the initialization lines
		seq_char = translation_table["".join([seq[i-2], seq[i-1], seq[i]])]
		i = i + 2

		for j in range(len(target)):  # columns
			target_char = target[j]
			j = j + 1  # skipping the initialization column
			# print(i, j)
			if seq_char == target_char:
				score = match
			else:
				score = mismatch

			dp_table[i][j] = max(
				# same frame
				dp_table[i-3][j-1] + score,
				dp_table[i][j-1] + gap,
				dp_table[i-3][j] + gap,

				# other frames
				dp_table[i-4][j-1] + fs,
				dp_table[i-2][j-1] + fs,
				0
				)

			# updating the max score and max score location
			if dp_table[i][j] > max_score:
				max_score = dp_table[i][j]
				max_score_loc = [(i, j)]
			elif dp_table[i][j] == max_score:
				max_score_loc.append((i,j))


	# print(max_score, max_score_loc)
	# print_dp(dp_table)
	# traceback from the locations of max score and find out where I went
	for start_i, start_j in max_score_loc:
		out_seq = ""
		out_target = ""
		diag_type = 3  # to tell me how many letters to add to out_seq when jumping

		while start_i > 1 or start_j > 0 or dp_table[start_i][start_j] != 0:
			# pdb.set_trace()
			i = start_i - 2  # because of the 2 initialization lines
			j = start_j - 1  # because of the 1 initialization column

			seq_char = translation_table["".join([seq[i-2], seq[i-1], seq[i]])]
			target_char = target[j]

			if seq_char == target_char:
				score = match
			else:
				score = mismatch

			######## within frame operations
			if dp_table[start_i][start_j] == dp_table[start_i-3][start_j-1] + score:  # match
				# I need to know what was my previous jump to know which letters to add to output
				if diag_type == 3:
					out_seq += seq[i] + seq[i-1] + seq[i-2]
				elif diag_type == 4:
					out_seq += seq[i] + seq[i-1] + seq[i-2] + seq[i-3]
				else:
					out_seq += seq[i] + seq[i-1]

				out_target += target_char
				diag_type = 3
				start_i, start_j = start_i-3, start_j-1


			elif dp_table[start_i][start_j] == dp_table[start_i][start_j-1] + score:  # deletion in seq
				out_seq += "---"
				out_target += target_char
				start_i, start_j = start_i, start_j-1


			elif dp_table[start_i][start_j] == dp_table[start_i-3][start_j] + score:  # insertion in seq
				out_seq += seq[i-2] + seq[i-1] + seq[i]
				out_target += "-"
				start_i, start_j = start_i-3, start_j


			####### frameshift operations
			elif dp_table[start_i][start_j] == dp_table[start_i-4][start_j-1] + fs:  # insertion shift, -1 frame
				if diag_type == 3:
					out_seq += seq[i] + seq[i-1] + seq[i-2]
				elif diag_type == 4:
					out_seq += seq[i] + seq[i-1] + seq[i-2] + seq[i-3]
				else:
					out_seq += seq[i] + seq[i-1]

				diag_type = 4
				out_target += target_char
				start_i, start_j = start_i-4, start_j-1


			elif dp_table[start_i][start_j] == dp_table[start_i-2][start_j-1] + fs:  # deletion shift, +1 frame
				if diag_type == 3:
					out_seq += seq[i] + seq[i-1] + seq[i-2]
				elif diag_type == 4:
					out_seq += seq[i] + seq[i-1] + seq[i-2] + seq[i-3]
				else:
					out_seq += seq[i] + seq[i-1]

				diag_type = 2
				out_target += target_char

				start_i, start_j = start_i-2, start_j-1

		print(out_seq[::-1])
		print(out_target[::-1])

	return dp_table


if __name__ == "__main__":

	if len(sys.argv) < 3:
		print("You need to give the query DNA seq and the target AA seq")
	seq = sys.argv[1]
	target = sys.argv[2]
	# print_dp(fs_aware_sw_codon_start(seq, target))
	table = fs_aware_sw_codon_end(seq, target)