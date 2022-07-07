import os
import sys
import random
import pdb


def read_fasta(fasta_file):
	sequences = dict()
	seqs = []
	seq_name = ""
	with open(fasta_file, "r") as infile:
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


in_dir = sys.argv[1]
files_percentage = int(sys.argv[2])/100
seq_percentage = int(sys.argv[3])/100

all_files = []
for f in os.listdir(in_dir):
	if f.endswith("fasta") or f.endswith("fa"):
		all_files.append(os.path.join(in_dir, f))

chosen = set()
while len(chosen) <= int(len(all_files)*files_percentage):
	chosen.add(random.choice(all_files))

# print("the chosen ones", len(chosen))
sequences = dict()
for f in chosen:
	seqs = read_fasta(f)
	# print("sequences: ", len(seqs))
	if len(seqs) == 1:
		for key, seq in seqs.items():
			sequences[f.split(os.sep)[1].split(".")[0]] = seq
	else:
		counter = 0
		selected_seqs = set()
		n_seqs = int(len(seqs) * seq_percentage)
		# print(n_seqs)
		if n_seqs == 0:
			n_seqs = 1
		keys = list(seqs.keys())
		while len(selected_seqs) < n_seqs:
			selected_seqs.add(random.choice(keys))
		# print("selected seqs:", len(selected_seqs))
		for key in selected_seqs:
			# the sequence will be named with the file name so I can match it to the graph when it aligns
			# pdb.set_trace()
			sequences[f.split(os.sep)[-1].split(".")[0] + "_" + key + "_" + str(counter)] = seqs[key]
			counter += 1

for key, seq in sequences.items():
	print(">" + key)
	print(seq)
