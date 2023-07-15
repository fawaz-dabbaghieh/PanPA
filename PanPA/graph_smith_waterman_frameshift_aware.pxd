from PanPA.Graph cimport Graph
from libcpp.vector cimport vector
from libcpp.string cimport string


cdef vector[string] align_to_graph_sw_fsa(Graph graph, str read, str read_name, bint print_dp,
                                          vector[int] sub_matrix, int gap_score, int fs_score,
                                          float min_id_score,
                                          vector[vector[vector[int]]] codon_translate_0base,
                                          vector[int] neucleotide_to_int) except *
