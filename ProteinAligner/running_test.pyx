# distutils: language=c++

import time
from ProteinAligner.Graph cimport Graph
from ProteinAligner.constants import all_linear_sub_matrices
from libcpp.vector cimport vector
import logging
import sys
import os
import gzip


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


def align_to_graph(graph_gfa, fasta_file, print_dp, gap_score, sub_matrix_name, output_file):
    cdef vector[int] sub_matrix
    cdef Graph graph = Graph(graph_gfa, paths=True)
    cdef int counter = 0
    # sub_matrix_name = "blosum62"
    cdef int i

    if not sub_matrix_name in all_linear_sub_matrices:
        print("Error! Check log file")
        logging.error("The substitution matrix {} is not present in the constants file".format(sub_matrix_name))
        sys.exit(1)
    else:
        for i in all_linear_sub_matrices[sub_matrix_name]:
            sub_matrix.push_back(i)

    # cdef Graph test_graph = Graph(gfa_file, paths=True)

    start = time.time()
    # print("starting to align {}".format(read))
    for read_name, read in read_fasta_gen(fasta_file):
        counter += 1
        read = "".join([x for x in read if x != "-"])
        min_id_score = 0.7
        graph.graph_align(read, read_name, print_dp, sub_matrix, gap_score, output_file, min_id_score)
    # print("It took {} for alignment".format(time.time() - start))
    # print("The graph size is {} mb".format(test_graph.__sizeof__()))
    print("{} alignments processed and it took {}".format(counter, time.time() - start))

def test_alignment(gfa_file, read, read_name, print_dp, gap_score, sub_matrix_name, output_file):
    cdef vector[int] sub_matrix
    # sub_matrix_name = "blosum62"
    cdef int i

    if not sub_matrix_name in all_linear_sub_matrices:
        print("Error! Check log file")
        logging.error("The substitution matrix {} is not present in the constants file".format(sub_matrix_name))
        sys.exit(1)
    else:
        for i in all_linear_sub_matrices[sub_matrix_name]:
            sub_matrix.push_back(i)

    cdef Graph test_graph = Graph(gfa_file, paths=True)

    start = time.time()
    # print("starting to align {}".format(read))
    min_id_score = 0.7

    test_graph.graph_align(read, read_name, print_dp, sub_matrix, gap_score, output_file, min_id_score)
    # print("It took {} for alignment".format(time.time() - start))
    # print("The graph size is {} mb".format(test_graph.__sizeof__()))


def test_alignment_ed(gfa_file, read, print_dp):
    cdef Graph test_graph = Graph(gfa_file, paths=True)
    start = time.time()
    print("starting to align {}".format(read))
    test_graph.graph_align_ed(read, print_dp)
    print("It took {} for alignment".format(time.time() - start))


def check_path(gfa_file, path):
    cdef Graph test_graph = Graph(gfa_file, paths=True)
    # print(path, test_graph.sorted)
    return(test_graph.path_seq(path))
