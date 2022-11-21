# distutils: language=c++

from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free
import random
import timeit
import time
from PanPA.constants import translation_table, all_linear_sub_matrices


cdef void print_dp_table(int all_seq_len, int read_len, int *dp_table):
    line = []
    cdef int i, row, j

    for i in range(read_len + 1):

        row = i * (all_seq_len + 1)

        for j in range(all_seq_len + 1):
            line.append(dp_table[row + j])
        print(line)
        line = []


def sw_def_cython(str read, vector[int] seq, bint print_table, sub_matrix_name = 'blosum62'):
    cdef vector[int] sub_matrix
    cdef int score
    for score in all_linear_sub_matrices[sub_matrix_name]:
        sub_matrix.push_back(score)

    cdef vector[int] read_int
    for c in read:
        read_int.push_back(ord(c))

    # cdef string sequence = bytes(sequence_py, 'utf-8')
    cdef int read_len = read_int.size()
    cdef int seq_len = seq.size()
    cdef int dimension = (read_len + 1) * (seq_len + 1)
    cdef int counter = 0
    # cdef vector[int] dp_table
    # dp_table.reserve(dimension)

    cdef int *dp_table = <int *> malloc(dimension * sizeof(int))

    cdef int gap_score = -3

    # cdef int match_score = 1

    cdef int max_score = 0
    # cdef int max_score_coordinates = 0
    # cdef vector[int] max_score
    cdef vector[int] max_score_coordinates
    cdef int row, left_cell, current_cell, above_cell, diagonal_cell, match, deletion, insertion, maximum, i, j, coord

    # cdef str out_read, out_seq

    # initializing vector with zeros, otherwise random values
    for i in range(dimension):
        dp_table[i] = 0

    # with nogil:
    for i in range(read_len):
        i += 1
        # because the dp table here is one dimension
        # I need to know in which row I am and every seq_len + 1 we reach a new row
        # i + 1 because I want to keep the first row 0, so there's an offset of one

        row = i * (seq_len + 1)
        for j in range(seq_len):
            counter += 1

            maximum = 0
            j += 1

            current_cell = row + j
            left_cell = current_cell - 1
            above_cell = current_cell - seq_len - 1
            diagonal_cell = above_cell - 1

            #"""
            match = dp_table[diagonal_cell] + sub_matrix[(read_int[i-1] - 65) * 27 + (seq[j-1] - 65)]

            # if read[i - 1] == sequence[j - 1]:
            #     match = dp_table[diagonal_cell] + match_score
            # else:
            #     match = dp_table[diagonal_cell] - match_score
            #"""
            #score = BLOSUM62[(read[i - 1],sequence[j - 1])]
            #match = dp_table[diagonal_cell] + score

            # match = dp_table[diagonal_cell] + (match_score if read[i] == sequence[j] else - match_score)
            deletion = dp_table[left_cell] + gap_score
            insertion = dp_table[above_cell] + gap_score

            maximum = max(match, deletion, insertion, 0)

            if max_score < maximum:
                max_score = maximum
                max_score_coordinates.clear()
                max_score_coordinates.push_back(current_cell)
            elif max_score == maximum:
                max_score_coordinates.push_back(current_cell)

            dp_table[current_cell] = maximum

# # traceback
#
#     for coord in max_score_coordinates:
#         max_score = dp_table[coord]
#         out_read = ""
#         out_seq = ""
#
#         # converting 1d to 2d coordinates
#         i = coord / (seq_len + 1)
#         j = coord % (seq_len + 1)
#         while ((i != 0) and (j != 0)) and (dp_table[coord] != 0):
#
#             current_cell = coord
#             left_cell = current_cell - 1
#             above_cell = current_cell - seq_len - 1
#             diagonal_cell = above_cell - 1
#
#             if (max_score == dp_table[diagonal_cell] + match_score) or (max_score == dp_table[diagonal_cell] - match_score):  # match or mismatch
#                 max_score = dp_table[diagonal_cell]
#                 coord = diagonal_cell
#                 i = coord / (seq_len + 1)
#                 j = coord % (seq_len + 1)
#                 out_read += read_py[i]
#                 out_seq += sequence_py[j]
#
#             elif max_score == dp_table[left_cell] + gap_score:  # insertion
#                 max_score = dp_table[left_cell]
#                 coord = left_cell
#                 i = coord / (seq_len + 1)
#                 j = coord % (seq_len + 1)
#                 out_read += "-"
#                 out_seq += sequence_py[j]
#
#
#             elif max_score == dp_table[above_cell] + gap_score:  # deletion
#                 max_score = dp_table[above_cell]
#                 coord = above_cell
#                 i = coord / (seq_len + 1)
#                 j = coord % (seq_len + 1)
#                 out_read += read_py[i]
#                 out_seq += "-"
#
#
#         # print(out_read[::-1])
#         # print(out_seq[::-1])
    print(counter)
    if print_table:
        print_dp_table(seq_len, read_len, dp_table)
    free(dp_table)

def run_sw_cython(s_len, r_len):
    amino_acids = list(set(translation_table.values()))
    amino_acids.remove("_")
    # print(amino_acids)
    seq = ''.join(random.choices(amino_acids, k=s_len))
    read = ''.join(random.choices(amino_acids, k=r_len))

    cdef vector[int] sequence
    for c in seq:
        sequence.push_back(ord(c))

    # t = timeit.timeit(lambda: sw_def_cython(read, sequence, 'blosum62'), number =20)
    # print("################", t/20)

    start = time.time()
    sw_def_cython(read, sequence, False, 'blosum62')
    print(time.time() - start)


def run_sw_cython_strings(read, seq):
    # amino_acids = list(set(translation_table.values()))
    # amino_acids.remove("_")
    # print(amino_acids)
    # seq = ''.join(random.choices(amino_acids, k=s_len))
    # read = ''.join(random.choices(amino_acids, k=r_len))

    cdef vector[int] sequence
    for c in seq:
        sequence.push_back(ord(c))

    # t = timeit.timeit(lambda: sw_def_cython(read, sequence, 'blosum62'), number =20)
    # print("################", t/20)

    start = time.time()
    sw_def_cython(read, sequence, True, 'blosum62')
    print(time.time() - start)