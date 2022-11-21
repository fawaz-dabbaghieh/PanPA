# distutils: language=c++
from libcpp cimport bool

cdef extern from "<algorithm>" namespace "std":
    Iter find[Iter, T](Iter first, Iter last, const T& value) except +

from libcpp.vector cimport vector

cdef class Node:
    # todo I need to change the identifier to a string, but this might have strong effects :(
    #    but maybe this would work if I remove the vector stuff and just use a list, I don't think using a vector
    #    here is faster anyway
    #    or keep the vector stuff but identifier as cpp string
    #    I think I will do this in the next release, after I write and submit the paper
    cdef int identifier
    cdef int coverage
    cdef str seq
    cdef vector[int] out_nodes
    cdef vector[int] in_nodes
    cdef set colors

    cdef void add_child(self, Node child) except *
    cdef void remove_child(self, Node child) except *
    cdef void add_parent(self, Node parent) except *
    cdef void remove_parent(self, Node parent) except *
    cdef bint is_child_of(self, Node parent) except *
    cdef bint is_parent_of(self, Node child) except *
