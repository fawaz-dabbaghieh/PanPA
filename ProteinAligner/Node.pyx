# distutils: language=c++
import sys
from libcpp.vector cimport vector
from libcpp cimport bool

# I decided to switch back to str and not use strings
# it introduced too much complexity and not much speed gain

cdef class Node:
    """
    A node object to store the a Node's information
    """

    def __init__(self, ident, seq=""):
        self.identifier = ident
        self.coverage = 1
        self.seq = seq
        self.out_nodes = set()
        self.in_nodes = set()
        self.colors = set()

    def __key(self):
        return self.identifier

    # read somewhere in cython documentation that special methods of a class need to bed def and not cdef
    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __sizeof__(self):
        total_size = 0
        total_size += sys.getsizeof(self.identifier)
        total_size += sys.getsizeof(self.seq)
        total_size += sys.getsizeof(self.out_nodes)
        total_size += sys.getsizeof(self.in_nodes)
        total_size += sys.getsizeof(self.colors)
        return total_size

    """
    Note: cdef functions that might have an assertion in them
    or potential error that need to be risen can't have void return
    Because otherwise the error will be created, reported, destroyed and
    the function will return from the exception point and the return statement will never be executed
    https://notes-on-cython.readthedocs.io/en/latest/misc.html
    Otherwise, if except * is added, then they can have void return and still report the error
    """
    cdef void add_child(self, Node child) except *:
        # assert isinstance(child, Node), "Children should be Node objects"
        # self.out_nodes.add(child.identifier)
        # child.in_nodes.add(self.identifier)
        # cdef vector[int].iterator found
        # found = find(self.out_nodes.begin(), self.out_nodes.end(), child.identifier)
        # if found == self.out_nodes.end():
        #     self.out_nodes.push_back(child.identifier)
        #     child.in_nodes.push_back(self.identifier)
        self.out_nodes.push_back(child.identifier)
        child.in_nodes.push_back(self.identifier)

    cdef void remove_child(self, Node child) except *:
        # assert child.identifier in self.out_nodes, "Child node {} is not one of " \
        #                                            "node's {} children".format(child.identifier, self.identifier)
        # self.out_nodes.remove(child.identifier)
        # child.in_nodes.remove(self.identifier)

        cdef vector[int].iterator found
        found = find(self.out_nodes.begin(), self.out_nodes.end(), child.identifier)

        assert found != self.out_nodes.end(), "Child node {} is not one of " \
                                                   "node's {} children".format(child.identifier, self.identifier)

        if found != self.out_nodes.end():
            self.out_nodes.erase(found)

        found = find(child.in_nodes.begin(), child.in_nodes.end(), self.identifier)
        if found != child.in_nodes.end():
            child.in_nodes.erase(found)

    cdef void add_parent(self, Node parent) except *:
        # self.in_nodes.add(parent.identifier)
        # parent.out_nodes.add(self.identifier)
        # cdef vector[int].iterator found
        # found = find(self.in_nodes.begin(), self.in_nodes.end(), parent.identifier)
        # if found == self.in_nodes.end():
        #     self.in_nodes.push_back(parent.identifier)
        #     parent.out_nodes.push_back(self.identifier)
        self.in_nodes.push_back(parent.identifier)
        parent.out_nodes.push_back(self.identifier)


    cdef void remove_parent(self, Node parent) except *:
        # assert parent.identifier in self.in_nodes, "Parent node {} is not one of " \
        #                                            "node's {} parents".format(parent.identifier, self.identifier)
        # self.in_nodes.remove(parent.identifier)
        # parent.out_nodes.remove(self.identifier)

        cdef vector[int].iterator found

        found = find(self.in_nodes.begin(), self.in_nodes.end(), parent.identifier)

        assert found != self.in_nodes.end(), "Parent node {} is not one of " \
                                                   "node's {} parents".format(parent.identifier, self.identifier)

        if found != self.in_nodes.end():
            self.in_nodes.erase(found)

        found = find(parent.out_nodes.begin(), parent.out_nodes.end(), self.identifier)
        if found != parent.out_nodes.end():
            parent.out_nodes.erase(found)

    cdef bint is_child_of(self, Node parent) except *:
        # if parent.identifier in self.in_nodes:
        #     return True
        # return False

        cdef vector[int].iterator found
        found = find(self.in_nodes.begin(), self.in_nodes.end(), parent.identifier)
        if found != self.in_nodes.end():
            return True
        return False

    cdef bint is_parent_of(self, Node child) except *:
        # if child.identifier in self.out_nodes:
        #     return True
        # return False

        cdef vector[int].iterator found
        found = find(self.out_nodes.begin(), self.out_nodes.end(), child.identifier)
        if found != self.out_nodes.end():
            return True
        return False
