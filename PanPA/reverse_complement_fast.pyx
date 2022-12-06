from libc.stdlib cimport free, malloc

cdef extern from "Python.h":
    const char * PyUnicode_AsUTF8AndSize(object, Py_ssize_t *size)

"""
This part done by Sanjay Srikakulam
"""

cpdef str reverse_complement(str seq):
    """TCreates reverse complement to a nucleotide sequence
    Args:
      seq: Any length nucleotide sequence
    Returns: Reverse complemented sequence
    """
    cdef char * basemap = [b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                           b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                           b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                           b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                           b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                           b'T', b'\0', b'G', b'\0', b'\0', b'\0', b'C', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'N',
                           b'\0', b'\0', b'\0', b'\0', b'\0', b'A', b'A', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0',
                           b'\0', b'\0', b'\0', b'\0', b'\0', b't', b'\0', b'g', b'\0', b'\0', b'\0', b'c', b'\0',
                           b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'\0', b'a', b'a']

    cdef Py_ssize_t i, seq_len
    cdef const char *seq_cstr = PyUnicode_AsUTF8AndSize(seq, &seq_len)
    cdef char *seq_pointer = &seq_cstr[0]
    cdef char *rev_comp = <char *> malloc(seq_len + 1)
    cdef str result

    for i in range(seq_len):
        rev_comp[seq_len - i - 1] = basemap[<int> seq_pointer[i]]

    result = rev_comp[:seq_len].decode('UTF-8')
    free(rev_comp)
    return result
