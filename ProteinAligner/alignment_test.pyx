# distutils: language=c++


import os
from ProteinAligner.Graph cimport Graph
from ProteinAligner.constants import all_linear_sub_matrices
from libcpp.vector cimport vector
from libcpp.string cimport string
from ProteinAligner.graph_smith_waterman cimport align_to_graph_sw


class AlignmentTest:
    def __init__(self):
        self.gfa = """S	1	MS	LN:i:2
L	1	+	4	+	0M
S	2	M	LN:i:1
L	2	+	4	+	0M
S	4	E	LN:i:1
L	4	+	5	+	0M
L	4	+	6	+	0M
S	5	PTPE	LN:i:4
L	5	+	14	+	0M
S	6	T	LN:i:1
L	6	+	8	+	0M
L	6	+	12	+	0M
S	8	QST	LN:i:3
L	8	+	14	+	0M
S	12	MA	LN:i:2
S	14	Q	LN:i:1
P	seq3	1+,4+,6+,8+,14+	0M,0M,0M,0M
P	seq1	2+,4+,5+,14+	0M,0M,0M
P	seq2	6+,12+	0M"""

    def test_alignment(self):
        cdef Graph graph
        cdef vector[string] alignments
        cdef vector[int] sub_matrix
        cdef int i
        cdef str seq, seq_name
        cdef bint print_db
        cdef int gap_score = -3

        for i in all_linear_sub_matrices["blosum62"]:
            sub_matrix.push_back(i)

        with open("testestest.gfa", "w") as outfile:
            outfile.write(self.gfa)
        graph = Graph("testestest.gfa")
        reads = {"seq1":"MEPTPEQ", "seq2":"MSTQSTQ"}
        seq = reads['seq1']
        seq_name = 'seq1'
        print_db = False
        alignments = align_to_graph_sw(graph, seq, seq_name, print_db, sub_matrix, gap_score, 0.5)
        assert alignments.size() == 1

        alignment = alignments[0].decode()
        alignment = alignment.split("\t")
        assert len(alignment) == 16
        assert alignment[0] == "seq1"
        assert alignment[5] == ">2>4>5>14"
        assert alignment[-1].split(":")[-1] == "7M"

        seq = reads['seq2']
        seq_name = 'seq2'
        print_db = False
        alignments = align_to_graph_sw(graph, seq, seq_name, print_db, sub_matrix, gap_score, 0.5)
        assert alignments.size() == 1
        alignment = alignments[0].decode()
        alignment = alignment.split("\t")
        assert len(alignment) == 16
        assert alignment[0] == "seq2"
        assert alignment[5] == ">1>4>6>8>14"
        assert alignment[-1].split(":")[-1] == "2M1D5M"

        os.remove("testestest.gfa")
        return True
