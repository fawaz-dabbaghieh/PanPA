from PanPA.node_test import NodeTest
'''
Wrapper functions to call the NodeTest class. Reason it was done like this because
The Node class modules are cdef and cannot be called directly from python
However, they can be called from another cython class or function that is def

NodeTest is a class and not separate functions, it's because pytest was complaining when 
cython def functions were called, but using a class it did not complain anymore.
Not sure why it gave this behavior but this was an easy workaround, I put all my tests in a class
and easily call it from these test_ python scripts
'''


def test_node_construction():
    tn = NodeTest()
    tn.test_node_construction()


def test_node_edges():
    tn = NodeTest()
    tn.test_node_edges()


def test_node_equality():
    tn = NodeTest()
    tn.test_equality()


def test_node_raising():
    tn = NodeTest()
    tn.test_error_rasing()
