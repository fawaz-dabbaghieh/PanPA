# distutils: language=c++

from PanPA.Node cimport Node
import pytest


# def cytest(func):
#     """
#     Wraps `func` in a plain Python function.
#     """
#     @functools.wraps(func)
#     def wrapped(*args, **kwargs):
#         bound = inspect.signature(func).bind(*args, **kwargs)
#         return func(*bound.args, **bound.kwargs)
#
#     return wrapped

# # @cytest
# def test_node_construction():
#     cdef Node new_node = Node(1, "what")
#     assert new_node.seq == b"what"
#     assert new_node.id == 1
#
# # @cytest
# def test_node_edges():
#     cdef Node new_node = Node(2, "what")
#     cdef Node child_node = Node(3, "?")
#     cdef Node parent_node = Node(1, "say")
#     new_node.add_child(child_node)
#     new_node.add_parent(parent_node)
#
#     assert child_node.id in new_node.out_nodes
#     assert new_node.id in parent_node.out_nodes
#     assert parent_node.id in new_node.in_nodes
'''
using the previous functions didn't work, because when they are imported into python
they are imported as builtin_function_or_method and pytest is not recognizing those as functions
So putting the test functions into a class and calling the class form test_*.py
works because pytest doesn't have this problem of not recognizing the function

However, this is a pytest problem, because running the test script with python interpreter 
gives normal behavior.
'''


class NodeTest:

    def test_node_construction(self):
        cdef Node new_node = Node(1, "ACGT")
        assert new_node.seq == "ACGT"
        assert new_node.identifier == 1

    def test_node_edges(self):
        cdef Node new_node = Node(2, "ACGT")
        cdef Node child_node = Node(3, "CGTA")
        cdef Node parent_node = Node(1, "TACG")
        new_node.add_child(child_node)
        new_node.add_parent(parent_node)

        assert new_node.is_parent_of(child_node)
        assert child_node.is_child_of(new_node)
        assert new_node.is_child_of(parent_node)
        assert parent_node.is_parent_of(new_node)

        new_node.remove_child(child_node)
        assert not new_node.is_parent_of(child_node)
        assert not child_node.is_child_of(new_node)

        new_node.remove_parent(parent_node)
        assert not new_node.is_child_of(parent_node)
        assert not parent_node.is_parent_of(new_node)

    def test_equality(self):
        cdef Node node1 = Node(1, "AAC")
        cdef Node node2 = Node(1, "AGA")
        cdef Node node3 = Node(2, "AAC")

        '''
        comparison is done on id
        the reason the test is based on id is because the way I construct my graphs from an MSA
        many nodes will have a single letter, but because the construction is done iteratively, new nodes will
        get unique ids
        '''
        assert node1 == node2
        assert node3 != node1

    def test_error_rasing(self):
        cdef Node node1 = Node(1, "GTC")
        # adding nodes function in the node class checks for the type of object given
        # if it's not a node object it will raise an error
        with pytest.raises(TypeError):
            node1.add_child(1)

        with pytest.raises(TypeError):
            node1.add_parent(2)

        cdef Node node2 = Node(2, "GGC")

        # node2 is not a parent or child of node 1, so there will be an assertion error
        # as I assert whether this is true or not in the remove_ functions of the node class
        with pytest.raises(AssertionError):
            node1.remove_parent(node2)

        with pytest.raises(AssertionError):
            node1.remove_child(node2)