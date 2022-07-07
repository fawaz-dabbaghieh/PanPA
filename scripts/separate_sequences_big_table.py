import gzip
import os
import logging
import sys
import pickle


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

if len(sys.argv) < 4:
    print("You need to give the big table pickled file and the sequences fasta file and file naming suffix")
    sys.exit()

in_table = sys.argv[1]
in_fasta = sys.argv[2]
suffix = sys.argv[3]


with open(in_table, "rb") as infile:
    big_table = pickle.load(infile)

all_seqs = dict()
for seq_name, seq in read_fasta_gen(in_fasta):
    all_seqs[seq_name] = seq


not_aligned = []
only_bwa = []
only_ga = []
only_panpa = []
only_bwa_ga = []
only_bwa_panpa = []
only_ga_panpa = []
all_three = []

for seq, info in big_table.items():
    if (info[A_BWA], info[A_GA], info[A_PANPA]) == (0,0,0):  # no alignments
        not_aligned.append(seq)
    elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (1,0,0):  # only BWA
        only_bwa.append(seq)
    elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (0,1,0):  # only GA
        only_ga.append(seq)
    elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (0,0,1):  # only PanPA
        only_panpa.append(seq)
    elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (1,1,0):  # BWA - GA
        only_bwa_ga.append(seq)
    elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (1,0,1):  # BWA - PanPA
        only_bwa_panpa.append(seq)
    elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (0,1,1):  # BA - PanPA
        only_ga_panpa.append(seq)
    elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (1,1,1):  # all of them
        all_three.append(seq)


with open("not_aligned" + "_" + suffix + ".fasta", "w") as out_file:
    for seq_name in not_aligned:
        out_file.write(f"{seq_name}\n{all_seqs[seq_name]}\n")

with open("only_bwa" + "_" + suffix + ".fasta", "w") as out_file:
    for seq_name in only_bwa:
        out_file.write(f"{seq_name}\n{all_seqs[seq_name]}\n")

with open("only_ga" + "_" + suffix + ".fasta", "w") as out_file:
    for seq_name in only_ga:
        out_file.write(f"{seq_name}\n{all_seqs[seq_name]}\n")

with open("only_panpa" + "_" + suffix + ".fasta", "w") as out_file:
    for seq_name in only_panpa:
        out_file.write(f"{seq_name}\n{all_seqs[seq_name]}\n")

with open("only_bwa_ga" + "_" + suffix + ".fasta", "w") as out_file:
    for seq_name in only_bwa_ga:
        out_file.write(f"{seq_name}\n{all_seqs[seq_name]}\n")

with open("only_bwa_panpa" + "_" + suffix + ".fasta", "w") as out_file:
    for seq_name in only_bwa_panpa:
        out_file.write(f"{seq_name}\n{all_seqs[seq_name]}\n")

with open("only_ga_panpa" + "_" + suffix + ".fasta", "w") as out_file:
    for seq_name in only_ga_panpa:
        out_file.write(f"{seq_name}\n{all_seqs[seq_name]}\n")

with open("all_three" + "_" + suffix + ".fasta", "w") as out_file:
    for seq_name in all_three:
        out_file.write(f"{seq_name}\n{all_seqs[seq_name]}\n")
