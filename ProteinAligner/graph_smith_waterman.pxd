from ProteinAligner.Graph cimport Graph
from libcpp.vector cimport vector
from libcpp.string cimport string


cdef vector[string] align_to_graph_sw(Graph graph, str read, str read_name, bint print_dp,
                             vector[int] sub_matrix, int gap_score, min_id_score) except *