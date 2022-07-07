import os
import sys
import matplotlib.pyplot as plt
from math import log

def read_fasta(fasta_file):
	sequences = dict()
	seqs = []
	seq_name = ""
	with open(os.path.join(in_dir, fasta_file), "r") as infile:
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
seq_files = []
n_seqs = []
n_samples = []
for f in os.listdir(in_dir):
	seqs = read_fasta(f)
	n_seqs.append(len(seqs))
	samples = set()
	for k in seqs.keys():
		name = "_".join(k.split("_")[0:2])  # the first two will have the accession number of the sample
		samples.add(name)
	n_samples.append(len(samples))
	seq_files.append(os.path.join(in_dir, f))


bins = int(1 + (3.22 * log(len(n_seqs))))
plt.hist(n_seqs, bins = bins)
plt.yscale('log')
plt.title("Sequence number distribution in clusters")
plt.xlabel("Num. of sequences")
plt.ylabel("Frequency")
plt.savefig("n_sequences_in_clusters.png", dpi=600)

plt.clf()

bins = int(1 + (3.22 * log(len(n_samples))))
plt.hist(n_samples, bins = bins)
plt.yscale('log')
plt.title("Assembly distribution in clusters")
plt.xlabel("Num. of assemblies")
plt.ylabel("Frequency")
plt.savefig("n_samples_in_clusters.png", dpi=600)
