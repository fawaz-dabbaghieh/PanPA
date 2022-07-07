import os
import sys
import gzip


def read_fasta(fasta_file_path):
    """
    read fasta file and return a dict of all sequences

    :param fasta_file_path: path to fasta file
    :return: a dictionary of seq_name:seq
    """

    sequences = dict()

    if not os.path.exists(fasta_file_path):
        print("file {} does not exist".format(fasta_file_path))
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


if len(sys.argv) < 2:
    print("You need to give a file with a list of MSA files")
    sys.exit()

in_dir = sys.argv[1]
list_of_files = []
with open(in_dir, "r") as in_file:
    for l in in_file:
        list_of_files.append(l.strip())

# log_file = open("warnings.log", "w")
for f in list_of_files:
    sequences = read_fasta(f)
    keys = list(sequences.keys())
    try:
        s_len = len(sequences[keys[0]])
    except IndexError:
        print(f"File {f} had an index error")
        continue
    for key, seq in sequences.items():
        try:
            assert len(seq) == s_len
        except AssertionError:
            print(f"File {f} is not an msa, sequences {key} does't seem to have the same size")
            break

print("All seems fine")
