# distutils: language=c++

from ProteinAligner.Graph cimport Graph
from constants import all_linear_sub_matrices
from libcpp.vector cimport vector

import time


class GraphTest:

    def test_graph_construction(self, gfa_file, paths):
        cdef Graph test_graph = Graph(gfa_file, paths=paths)
        # test_graph.read_gfa(gfa_file, paths=True)
        assert len(test_graph) == 15
        assert len(test_graph.paths) == 3

        assert "seq1" in test_graph.paths
        assert "seq2" in test_graph.paths
        assert "seq3" in test_graph.paths


    def test_graph_writing(self, gfa_file):
        cdef Graph test_graph = Graph(gfa_file, paths=True)
        # test_graph.read_gfa(gfa_file, paths=True)
        test_graph.write_gfa("writing_testing.gfa")


    def test_graph_compacting(self, gfa_file):
        cdef Graph test_graph = Graph(gfa_file, paths=True)
        # test_graph.read_gfa(gfa_file, paths=True)
        n_nodes = len(test_graph)

        test_graph.compact()
        print(test_graph.paths)
        test_graph.write_gfa("compacted_example.gfa")
        assert n_nodes > len(test_graph)


    def test_top_sorting(self, gfa_file):
        cdef Graph test_graph = Graph(gfa_file, paths=True)
        # test_graph.read_gfa(gfa_file, paths=True)
        assert len(test_graph.sorted) == len(test_graph)


