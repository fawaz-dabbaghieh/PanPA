# distutils: language=c++

import sys
import logging
from PanPA.Graph cimport Graph
from PanPA.Node cimport Node


cdef class Alignment:
    def __init__(self, read_name, read_len, alignment_score):
        self.read_name = read_name
        self.read_len = read_len
        self.alignment_score = alignment_score
        self.path = []
        self.info = []
        self.n_matches = 0
        self.n_mismatches = 0
        self.n_indels = 0
        self.id_score = 0
        self.gaf = ""


    cdef void prepare_aa_gaf(self, Graph graph) except *:
        """
        Outputs alignment in GAF format, either to std output or to a file
        :param graph: is a graph object
        :param output_file: is an opened file with "w" object
        :param stdout: boolean to whether to write to stdout or not
        """
        if not self.info:
            logging.error("No GAF for this alignment, the info is {}".format(self.read_name))
            sys.exit()

        cdef Node node

        # todo I can just loop backwards instead of reversing everytime
        self.info.reverse()
        self.path.reverse()

        # print(f"I am in prepare gaf and info is reverse and it is {self.info}")
        # should never reach this point
        # if (output_file is None) and (stdout is False):
        #     print("You need to either give an output alignment file or True for stdout to print the alignment")
        #     sys.exit()

        gaf_string = [self.read_name, self.read_len, self.info[0]["read_pos"], self.info[-1]["read_pos"] + 1, "+"]

        # adding the path
        path = []
        path_seq_len = 0
        for n in self.path:
            if not path:
                node = graph.nodes[n]
                # path_seq_len += len(node.seq)
                path.append(n)
            elif not path[-1] == n:
                # path_seq_len += len(node.seq)
                path.append(n)
            else:
                continue

        for p in path:
            node = graph.nodes[p]
            path_seq_len += len(node.seq)
        path = "".join([">" + str(x) for x in path])

        gaf_string.append(path)
        # for n in alignment.path:
        #     path_seq_len += len(self.graph.nodes[n].seq)
        gaf_string.append(path_seq_len)

        # for now I'll keep it related to the number of nodes in path
        gaf_string += [self.info[0]["node_pos"], self.info[0]["node_pos"] + len(self.path)]


        gaf_string.append(int(self.n_matches))
        # I think I need to add +1 here because it's [start, end[
        # so like when you say in python3 range(1,5) you get a count until 4
        gaf_string.append(len(self.info))
        gaf_string.append(255)

        # identity = round(self.n_matches/len(self.info), 5)
        # I added some tags, as for alignment score (the score in the DP table)
        # dv for per-base sequence divergence (I think it's 1 - identity)
        # id for sequence identity
        # Maybe add tp:A:(P, S) for primary or secondary alignment
        # cm:i number of minimizers (maybe see how many seeds from the read hit that graph)
        # NM:i total numbers of mismatches and gaps (indels)
        # print(f"calculating the id_score from n matches {self.n_matches}")
        self.id_score = self.n_matches / float(len(self.info))
        gaf_string.append(f"NM:i:{self.n_indels + self.n_mismatches}")
        gaf_string.append(f"AS:i:{str(self.alignment_score)}")
        gaf_string.append(f"dv:f:{str(round(1-self.id_score, 4))}")
        gaf_string.append(f"id:f:{str(round(self.id_score, 4))}")
        # previous = (alignment.info[0]["cigar"], 1)
        previous = []
        cigar = "cg:Z:"
        cigar_symbols = ["I", "D", "=", "X"]
        # pdb.set_trace()
        for item in self.info:
            if not previous:
                previous = [item["type"], 1]

            elif item["type"] == previous[0]:
                previous[1] += 1
            else:
                cigar += str(previous[1]) + cigar_symbols[previous[0]]  # first letter (M, I, D)
                previous = [item["type"], 1]
        cigar += str(previous[1]) + cigar_symbols[previous[0]]

        gaf_string.append(cigar)

        # adding graph name as an extra tag
        gaf_string.append(f"gid:Z:{graph.name}")

        self.gaf = "\t".join([str(x) for x in gaf_string])


    cdef void prepare_dna_gaf(self, Graph graph) except *:
        """
        Outputs alignment in GAF format, either to std output or to a file
        :param graph: is a graph object
        :param output_file: is an opened file with "w" object
        :param stdout: boolean to whether to write to stdout or not
        """
        if not self.info:
            logging.error("No GAF for this alignment, the info is {}".format(self.read_name))
            sys.exit()

        cdef Node node

        # todo I can just loop backwards instead of reversing everytime
        self.info.reverse()
        self.path.reverse()

        # print(f"I am in prepare gaf and info is reverse and it is {self.info}")
        # should never reach this point
        # if (output_file is None) and (stdout is False):
        #     print("You need to either give an output alignment file or True for stdout to print the alignment")
        #     sys.exit()

        gaf_string = [self.read_name, self.read_len, self.info[0]["read_pos"], self.info[-1]["read_pos"] + 1, "+"]

        # adding the path
        path = []
        for n in self.path:
            if not path:
                node = graph.nodes[n]
                # path_seq_len += len(node.seq)
                path.append(n)
            elif not path[-1] == n:
                # path_seq_len += len(node.seq)
                path.append(n)
            else:
                continue

        path_seq_len = 0
        for p in path:
            node = graph.nodes[p]
            path_seq_len += len(node.seq)
        path = "".join([">" + str(x) for x in path])

        gaf_string.append(path)
        # for n in alignment.path:
        #     path_seq_len += len(self.graph.nodes[n].seq)
        gaf_string.append(path_seq_len)

        # for now I'll keep it related to the number of nodes in path
        gaf_string += [self.info[0]["node_pos"], self.info[0]["node_pos"] + len(self.path)]

        gaf_string.append(int(self.n_matches))
        # I think I need to add +1 here because it's [start, end[
        # so like when you say in python3 range(1,5) you get a count until 4
        gaf_string.append(len(self.info))
        gaf_string.append(255)

        # identity = round(self.n_matches/len(self.info), 5)
        # I added some tags, as for alignment score (the score in the DP table)
        # dv for per-base sequence divergence (I think it's 1 - identity)
        # id for sequence identity
        # Maybe add tp:A:(P, S) for primary or secondary alignment
        # cm:i number of minimizers (maybe see how many seeds from the read hit that graph)
        # NM:i total numbers of mismatches and gaps (indels)
        # print(f"calculating the id_score from n matches {self.n_matches}")
        self.id_score = self.n_matches / float(len(self.info))
        gaf_string.append(f"NM:i:{self.n_indels + self.n_mismatches}")
        gaf_string.append(f"AS:i:{str(self.alignment_score)}")
        gaf_string.append(f"dv:f:{str(round(1-self.id_score, 4))}")
        gaf_string.append(f"id:f:{str(round(self.id_score, 4))}")
        # previous = (alignment.info[0]["cigar"], 1)
        previous = []
        cigar = "cg:Z:"
        cigar_symbols = ["I", "D", "=", "X"]
        # pdb.set_trace()
        for item in self.info:
            if not previous:
                previous = [item["type"], 1]

            elif item["type"] == previous[0]:
                previous[1] += 1
            else:
                cigar += str(previous[1]) + cigar_symbols[previous[0]]  # first letter (M, I, D)
                previous = [item["type"], 1]
        cigar += str(previous[1]) + cigar_symbols[previous[0]]

        gaf_string.append(cigar)

        # adding graph name as an extra tag
        gaf_string.append(f"gid:Z:{graph.name}")

        self.gaf = "\t".join([str(x) for x in gaf_string])
