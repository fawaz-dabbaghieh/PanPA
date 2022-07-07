import os
import sys
import gzip

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
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

reverse = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def read_fasta(fasta_file_path):
    """
    read fasta file and return a dict of all sequences

    :param fasta_file_path: path to fasta file
    :return: a dictionary of seq_name:seq
    """

    sequences = dict()

    if not os.path.exists(fasta_file_path):
        print("file does not exist")
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


def get_aa(seq, start, end, strand):
    if strand == "-":  # because reverse complement
        tmp_end = len(seq) - start + 1
        tmp_start = len(seq) - end
        start = tmp_start
        end = tmp_end

        seq_rev = "".join([reverse[c] for c in seq[::-1]])
        sub_seq = seq_rev[start:end]
    else:
        start = start - 1
        sub_seq = seq[start:end]
    
    assert len(sub_seq) % 3 == 0  # make sure it's the full protein

    aa = []
    for i in range(0, len(sub_seq), 3):
        aa.append(translation_table[sub_seq[i:i+3]])

    return "".join(aa[:-1])


if len(sys.argv) < 4:
    print("you need to give an assembly, info table, target protein patric ID")
    sys.exit()


assembly = sys.argv[1]
strain = assembly.split(os.sep)[-1].split(".fna")[0]
table = sys.argv[2]
target = sys.argv[3]
START = 9
END = 10
STRAND = 11
CONTIG = 2

all_targets = []
with open(table, "r") as infile:
    for l in infile:
        if target in l:
            all_targets.append(l.split("\t"))

sequences = read_fasta(assembly)

for t in all_targets:
    start = int(t[START])
    end = int(t[END])
    strand = t[STRAND]
    contig = t[CONTIG]
    for seq_id, seq in sequences.items():
        if contig in seq_id:  # found the contig where our target is
            print(f">{strain}_{contig}_{target}_{start}_{end}")
            print(get_aa(seq.upper(), start, end, strand))
