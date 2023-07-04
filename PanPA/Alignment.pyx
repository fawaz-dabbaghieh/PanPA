# distutils: language=c++

import sys
import json
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

    cdef void add_alignment(self, int node_id, int node_letter_int, int node_pos, int read_letter_int,
                            int read_pos) except *:
        """
        add each letter alignment to the info list which will be converted later to GAF format

        for read_letter_int and graph_seq_letter_int, if -1 then gap, otherwise I need to convert the int to string
        with chr(read_letter + 65)
        """
        if node_letter_int == -1:
            node_str = "_"
        else:
            node_str = chr(node_letter_int + 65)

        if read_letter_int == -1:
            read_str = "_"
        else:
            read_str = chr(read_letter_int + 65)

        if node_str == "_":
            cigar = "Insertion"
            self.n_indels += 1
        elif read_str == "_":
            cigar = "Deletion"
            self.n_indels += 1
        elif node_str == read_str:
            cigar = "Match"
            self.n_matches += 1
        elif node_str != read_str:
            cigar = "Mismatch"
            self.n_mismatches += 1

        else:
            print("one condition for cigar was not processed, exit")
            sys.exit()


        self.info.append({"node_id": node_id, "node_pos": node_pos, "node_str": node_str,
                          "read_pos": read_pos, "read_str": read_str, "cigar": cigar})

        self.path.append(node_id)


    cdef void prepare_gaf(self, Graph graph, min_id_score) except *:
        """
        Outputs alignment in GAF format, either to std output or to a file
        :param graph: is a graph object
        :param min_id_score: identity score between 0 and 1 to filter the alignment against
        """
        if not self.info:
            logging.error("No GAF for this alignment, the info is {}".format(self.read_name))
            sys.exit()

        cdef Node node

        # todo I can just loop backwards instead of reversing everytime
        self.info.reverse()
        self.path.reverse()

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
        cigar_symbols = {"Match":"=", "Mismatch":"X", "Deletion":"D", "Insertion":"I"}
        for item in self.info:
            if not previous:
                previous = [item["cigar"], 1]

            elif item["cigar"] == previous[0]:
                previous[1] += 1
            else:
                cigar += str(previous[1]) + cigar_symbols[previous[0]]  # first letter (M, I, D)
                previous = [item["cigar"], 1]
        cigar += str(previous[1]) + cigar_symbols[previous[0]]

        gaf_string.append(cigar)

        # adding graph name as an extra tag
        gaf_string.append(f"gid:Z:{graph.name}")

        self.gaf = "\t".join([str(x) for x in gaf_string])
        # return "\t".join([str(x) for x in gaf_string]) + "\n"
        # very hacky
        # if stdout:
        #     print("\t".join([str(x) for x in gaf_string]))
        # else:
        #     output_file.write("\t".join([str(x) for x in gaf_string]))
        #     output_file.write("\n")


    def print_alignment(self):
        print(json.dumps(self.info))