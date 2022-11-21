# distutils: language=c++
from libcpp.vector cimport vector


cdef class Graph:

    cdef dict nodes, paths
    # cdef dict column_pos
    # cdef vector[vector[int]] in_nodes_map
    cdef str name

    cdef vector[int] node_ends, intervals, in_nodes
    # self.node_ends, self.column_pos = dict()

    # maybe change lists to int * and then malloc later when initializing
    cdef list sorted
    cdef vector[int] j_node, j_pos
    cdef int seq_len
    # might remove for saving memory later
    # instead, when having a j, just access the node's string for the letter
    cdef vector[int] all_seqs

    cdef void prepare_graph(self) except *
    cdef void sort(self) except *
    cdef void add_paths(self) except *
    cdef void nodes_info(self, output_file) except *
    cdef bint merge_end(self, int n) except *
    cdef void compact(self) except *
    cdef void top_sorting(self) except *

    cdef void read_gfa(self, str gfa_path, paths=*) except *
    cdef void write_gfa(self, str gfa_path) except *
    cdef str path_seq(self, list path)
    # should probably add a remove node function
    # that removes a node safely (i.e. removes all the associated edges properly)
    # using remove parent and remove child functions from the node object
    # then removes the node from the node list

