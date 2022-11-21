from PanPA.Graph cimport Graph


cdef class Alignment:
    cdef list path, info
    cdef str read_name
    cdef int read_len
    cdef int alignment_score
    cdef int n_matches
    cdef int n_mismatches
    cdef int n_indels
    cdef float id_score
    cdef str gaf

    cdef void add_alignment(self, int node_id, int node_letter_int, int node_pos, int read_letter_int,
                            int read_pos) except *

    cdef void prepare_gaf(self, Graph graph, min_id_score) except *
