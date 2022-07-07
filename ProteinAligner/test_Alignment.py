from ProteinAligner.alignment_test import AlignmentTest


def test_alignment():
    test = AlignmentTest()
    assert test.test_alignment() == True
