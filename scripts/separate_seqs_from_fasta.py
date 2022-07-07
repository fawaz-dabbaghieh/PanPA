import gzip
import os
import logging
import sys


def read_fasta_gen(fasta_file_path):
    """
    A generator function that reads one read at a time
    Can be used for big FASTA files to not keep them in memory
    :param fasta_file_path: path to fasta file
    :yield: a tuple of sequence id and sequence
    """

    if not os.path.exists(fasta_file_path):
        logging.error("file {} does not exist".format(fasta_file_path))
        sys.exit()

    if fasta_file_path.endswith("gz"):
        fasta_file = gzip.open(fasta_file_path, "rt")
    else:
        fasta_file = open(fasta_file_path, "r")

    seqs = []
    seq_name = ""
    for line in fasta_file:
        line = line.strip()
        if not line:  # empty line
            continue

        if line.startswith(">"):
            if len(seqs) != 0:  # there was a sequence before
                yield seq_name, "".join(seqs)
                seq_name = line[1:]
                seqs = []
            else:
                seq_name = line[1:]
        else:
            seqs.append(line)

    # last sequence
    if seqs:
        yield seq_name, "".join(seqs)



if len(sys.argv) < 3:
    print("You need to give the fasta file and a text file with seq names")
    sys.exit()

fasta_file = sys.argv[1]
seqs_file = sys.argv[2]
to_output = set()
with open(seqs_file, "r") as infile:
    for l in infile:
        to_output.add(l.strip())

for seq_name, seq in read_fasta_gen(fasta_file):
   if seq_name in to_output:
        print(">" + seq_name)
        print(seq)
