from ProteinAligner.graph_test import GraphTest


def test_graph_construction(gfa_file="testing/example_seqs.gfa", paths=True):
    gt = GraphTest()
    gt.test_graph_construction(gfa_file, paths)


def test_graph_writing(gfa_file="testing/example_seqs.gfa"):
    gt = GraphTest()
    gt.test_graph_writing(gfa_file)


def test_graph_compacting(gfa_file="testing/example_seqs.gfa"):
    gt = GraphTest()
    gt.test_graph_compacting(gfa_file)


def test_top_sorting(gfa_file="testing/example_seqs.gfa"):
    gt = GraphTest()
    gt.test_top_sorting(gfa_file)
