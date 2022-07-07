from ProteinAligner.Graph cimport Graph
from libcpp.vector cimport vector
from libcpp.map cimport map

cpdef void align_to_graph_ed(Graph graph, str read, bint print_dp) except *