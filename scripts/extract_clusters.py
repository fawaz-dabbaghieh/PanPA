import gzip
import os
import sys
import logging


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
                # seq_name without >
                yield seq_name, "".join(seqs)
                seq_name = line[1:]
                seqs = []
            # because mmseqs2 output has the name of the cluster as a fasta sequence id
            # then no sequence after that, so there will be two id lines consecutively
            # and I want to return that name because it's the cluster id so i konw where my cluster starts
            elif seq_name:
                yield seq_name, ""
            else:    
                seq_name = line[1:]
        else:
            seqs.append(line)

    # last sequence
    if seqs:
        yield seq_name, "".join(seqs)


if len(sys.argv) < 4:
    print("You need to provide <list_of_cluster_names.txt> <mmseqs_fasta_output_with_all_sequences.fasta> <out_dir>")
    print("where the list is the name of each cluster id which is usually the sequence that represent that cluster and each on one line")
    print("the second one is the mmseqs2 output which has all the sequence and the cluster names")
    sys.exit()
list_of_cluster_file = sys.argv[1]
in_clusters_fasta = sys.argv[2]
out_dir = sys.argv[3]

list_of_clusters = set()
with open(list_of_cluster_file, "r") as in_file:
    for l in in_file:
        list_of_clusters.add(l.strip())

previous_cluster_name = ""
previous_cluster = dict()
for seq_name, seq in read_fasta_gen(in_clusters_fasta):

    if not seq:
        # we are reading a cluster name because mmseqs introduces the name of the cluster with an empty
        # fasta seq, so there's a > line with name of cluster then immidiately another > line with the first seq
        # in the cluster.
        # if we had a previous then we need to write it out because we are starting a new cluster
        if previous_cluster_name and previous_cluster_name in list_of_clusters:
            # we have a full cluster and we need to output it
            name = os.path.join(out_dir, previous_cluster_name + ".fasta")
            with open(name, "w") as out_file:
                for s_id, s in previous_cluster.items():
                    out_file.write(">" + s_id + "\n")
                    out_file.write(s + "\n")
                # I wrote the previous cluster in a file and emptying the cluster for the new one
            previous_cluster = dict()
            # print(previous_cluster_name)
            # del previous_cluster_name
            previous_cluster_name = seq_name  # new cluster name

        else:  # we don't have a previous cluster name, so it's the first cluster
            previous_cluster_name = seq_name
            previous_cluster = dict()
        # else:
        #     # a cluster I don't want, so I write the previous cluster if there was one
        #     if previous_cluster_name:
        #         name = os.path.join(out_dir, previous_cluster_name.split(":")[0] + ".fasta")
        #         with open(name, "w") as out_file:
        #             for s_id, s in previous_cluster.items():
        #                 out_file.write(">" + s_id + "\n")
        #                 out_file.write(s + "\n")
        #     previous_cluster = dict()
        #     previous_cluster_name = ""
            
    else:  # it's a normal sequence and I'll keep it for now until I've read the whole cluster
        previous_cluster[seq_name] = seq

if previous_cluster_name and previous_cluster_name in list_of_clusters:  # last cluster
    name = os.path.join(out_dir, previous_cluster_name.split(":")[0] + ".fasta")
    with open(name, "w") as out_file:
        for s_id, s in previous_cluster.items():
            out_file.write(">" + s_id + "\n")
            out_file.write(s + "\n")
