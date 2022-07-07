# distutils: language=c++

import logging
import sys
import time
from ProteinAligner.Graph cimport Graph
from ProteinAligner.Alignment cimport Alignment
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free


cdef void print_dp_table(int all_seq_len, int read_len, int *dp_table):
    line = []
    cdef int i, row, j

    for i in range(read_len + 1):

        row = i * (all_seq_len + 1)

        for j in range(all_seq_len + 1):
            line.append(dp_table[row + j])
        print(line)
        line = []
    print("\n")


cpdef void align_to_graph_ed(Graph graph, str read, bint print_dp) except *:
    """
    getting rid of nested vectors
    [0, 1, 2, 2, 4, 4, 3, 5]
    [0, 0, 0, 0, 0, 2, ]  (two entries for each node i, i+1 for the interval)
    
    [] one entry for each node, giving starting for each interval


    In this part, I am assigning pointers to structures in the graph object
    otherwise, if I access them with graph.object, cython is adding extra python interaction lines
    (yellow lines in the annotated html file) that can cause some slow down
    
    when assigning for example graph.j_pos to a vector I assigned, it's not being copied (graph.j_pos and j_pos get
    the same id) but then cython doesn't add these extra checks
    
    I think the checks were related to whether the graph object has that attribute or not
    """
    #########################################################

    # cdef vector[vector[int]] in_nodes_map
    # in_nodes_map = graph.in_nodes_map
    cdef vector[int] intervals
    intervals = graph.intervals
    cdef vector[int] in_nodes_vec
    in_nodes_vec = graph.in_nodes

    cdef vector[int] j_node, j_pos, node_ends
    node_ends = graph.node_ends
    j_node = graph.j_node
    j_pos = graph.j_pos

    cdef vector[int] in_nodes
    cdef size_t in_nodes_size

    cdef vector[int] all_seq_as_int = graph.all_seqs
    #########################################################

    #########################################################
    """
    Defining some variables I need for the algorithm
    """
    cdef int which_node = -1  # counter to tell me in which node I am in the column loop

    # keeping track of all possible max scores in the whole dp_table
    cdef int global_min = 50000  # some random big value
    cdef vector[int] global_min_coord  # cell coordinates of all max scores (if more than on optimal alignment)

    cdef int in_nodes_max, in_nodes_min_coord  # which in node gives the best score

    cdef int n, dimensions, i, j, left_cell, current_cell, above_cell, diagonal_cell, row, l
    cdef int match_miss, deletion, insertion, local_min, previous_node, previous_j
    # cdef int gap_score = -3
    cdef int current_node_id, current_node_end, current_node_pos

    cdef int counter = 0  # for testing
    #########################################################

    #########################################################
    """
    Converting the read to integers
    declaring the dp table and filling it with 0 values
    
    and filling traceback table
    """
    cdef vector[int] read_as_int
    cdef int *dp_table  # dp table is an array with memory allocated
    cdef int read_len, graph_seq_len

    read_len = len(read)
    graph_seq_len = graph.all_seqs.size()

    dimensions = (read_len + 1) * (graph_seq_len + 1)
    dp_table = <int *> malloc(dimensions * sizeof(int))
    # initializing the vector

    for i in range(dimensions):
        dp_table[i] = 0

    for c in read:
        read_as_int.push_back(ord(c) - 65)

    cdef int *traceback_table
    traceback_table = <int *> malloc(dimensions * sizeof(int))
    for i in range(dimensions):
        traceback_table[i] = 0
    #########################################################
    # with nogil:
    for i in range(read_len):
        # adding one to keep first row 0
        i += 1

        which_node = -1  # resetting the node counter

        # because the dp table is 1d, I need to convert between i,j to 1d coordinates
        row = i * (graph_seq_len + 1)
        # I want to align the whole sequence but to part of the graph
        # so the first row is all 0, but first column is incremented
        dp_table[row] = i - 1

        for j in range(graph_seq_len):
            # print("i {}, j {}".format(i, j))

            # counter += 1
            maximum = 0
            # adding 1 to keep first column 0
            j += 1

            """
            I have 3 situations:
            1- I am at the left border of a node and the node doesn't have in nodes
            2- I am not at the left border of the node
            For (1 and 2) I just need to take i-1 and j-1, so can be grouped together
            
            3- I am at the left border and the node has in nodes
            This one then need the column jumping, I look at the in nodes, and for each node I take the j from the
            node_end dictionary
            """
            in_nodes.clear()
            in_nodes_size = 0

            # current_node_id = j_node[j-1]  # which node corresponds to the current j
            current_node_pos = j_pos[j-1]  # where in that node am I
            if current_node_pos == 0:  # we are at a new node
                which_node += 1

            for l in range(intervals[which_node], intervals[which_node + 1]):
                in_nodes.push_back(in_nodes_vec[l])
            # in_nodes = in_nodes_map[which_node]
            in_nodes_size = in_nodes.size()

            # print(i, j, current_node_id, which_end, in_nodes)
            # todo maybe not use traceback matrix, because it's adding around 20% time on the base algorithm
            #   because now I need to fill a traceback matrix for each cell
            #   without it, I only need to do max() calculations while doing the traceback


            # 1 and 2 situation
            if (current_node_pos != 0) or (current_node_pos == 0 and in_nodes_size == 0):

                current_cell = row + j
                left_cell = current_cell - 1
                above_cell = current_cell - graph_seq_len - 1
                diagonal_cell = above_cell - 1

                # because sub_matrix is linearized over all the alphabet (26 characters)
                # both read_as_int and all_seq_as_int are integers with A being 0 and Z 25
                # so accessing a certain cell in the 2d substitution matrix corresponds to
                # first_letter * 26 + second_letter
                if read_as_int[i-1] == all_seq_as_int[j-1]:
                    match_miss = dp_table[diagonal_cell]
                else:
                    match_miss = dp_table[diagonal_cell] + 1

                deletion = dp_table[left_cell] + 1
                insertion = dp_table[above_cell] + 1

                local_min = min(match_miss, deletion, insertion)

                # filling the traceback table
                if local_min == match_miss:
                    traceback_table[current_cell] = diagonal_cell
                elif local_min == deletion:
                    traceback_table[current_cell] = left_cell
                elif local_min == insertion:
                    traceback_table[current_cell] = above_cell

                # # updating global max if needed
                # if local_min < global_min:
                #     global_min = local_min
                #     global_min_coord.clear()
                #     global_min_coord.push_back(current_cell)
                # elif local_min == global_min:
                #     global_min_coord.push_back(current_cell)


                dp_table[current_cell] = local_min
                # print(current_cell, left_cell, above_cell, diagonal_cell)
                # print(i, j, chr(read_as_int[i-1] + 65), chr(all_seq_as_int[j-1] + 65), match_miss, deletion, insertion, local_min)


            # situation 3
            elif current_node_pos == 0:

                # todo I don't think I need this if, because both zero is in first situation
                if in_nodes_size != 0:


                    in_nodes_min = 100
                    in_nodes_min_coord = 0
                    # current_cell, above_cell and insertion don't depend on the incoming edges
                    # so only need to be calculated once
                    current_cell = row + j
                    above_cell = current_cell - graph_seq_len - 1

                    insertion = dp_table[above_cell] + 1

                    for previous_node in in_nodes:

                        # node_ends tells me where to jump in the dp table
                        # so instead of j-1 I get some other previous_j
                        previous_j = node_ends[previous_node]

                        left_cell = row + previous_j
                        diagonal_cell = left_cell - graph_seq_len - 1

                        if read_as_int[i - 1] == all_seq_as_int[j-1]:

                            match_miss = dp_table[diagonal_cell]
                        else:
                            match_miss = dp_table[diagonal_cell] + 1
                        #
                        # print(i, j, read[i - 1], node.seq[current_node_pos],
                        #       sub_matrix[read_as_int[i-1] + all_seq_as_int[j-1]], in_nodes_size)

                        deletion = dp_table[left_cell] + 1

                        # updating in_nodes_min which will tell me which jump had the best score
                        local_min = min(match_miss, deletion, insertion)

                        if local_min <  in_nodes_min:
                            # filling the traceback table with the potential max
                            if local_min == match_miss:
                                in_nodes_min_coord = diagonal_cell
                            elif local_min == deletion:
                                in_nodes_min_coord = left_cell
                            elif local_min == insertion:
                                in_nodes_min_coord = above_cell

                            in_nodes_min = local_min

                        # print(i, j, current_cell, left_cell, above_cell, diagonal_cell, local_max)

                    # updating the global max and its coordinates
                    # if in_nodes_min < global_min:
                    #     global_min = local_min
                    #     global_min_coord.clear()
                    #     global_min_coord.push_back(current_cell)
                    # elif in_nodes_min == global_min:
                    #     global_min_coord.push_back(current_cell)

                    traceback_table[current_cell] = in_nodes_min_coord
                    dp_table[current_cell] = in_nodes_min

                    # print(current_cell, left_cell, above_cell, diagonal_cell)
                    # print(i, j, chr(read_as_int[i-1] + 65), chr(all_seq_as_int[j-1] + 65), match_miss, deletion, insertion, in_nodes_min)

    # print(counter)
    # getting the minimum in the last row
    last_row_begin = read_len * (graph_seq_len + 1)
    last_row_end = dimensions
    for i in range(last_row_begin, last_row_end):
        if dp_table[i] < global_min:
            global_min = dp_table[i]
            global_min_coord.clear()
            global_min_coord.push_back(i)
        elif dp_table[i] == global_min:
            global_min_coord.push_back(i)

    if print_dp:
        print("This is the dp table\n")
        print_dp_table(graph_seq_len, read_len, dp_table)
        print("This is the traceback table\n")
        print_dp_table(graph_seq_len, read_len, traceback_table)

    ###############################################################################################
    # traceback now
    # todo I need to give this function the read name as well, I can do this later
    cdef int coord, back_coord
    cdef Alignment alignment
    cdef int traceback_max, back_j, back_i

    print("what", global_min, global_min_coord)
    for coord in global_min_coord:
        # max_score = dp_table[coord]
        alignment = Alignment("seq1", read_len)

        # traceback_max = dp_table[coord]
        i = coord / (graph_seq_len + 1)
        j = coord % (graph_seq_len + 1)
        # i and j can be used as is to access the corresponding letter in

        while ((i != 0) and (j != 0)):
            # print("current: coord, i, j", coord, i, j)
            back_coord = traceback_table[coord]
            back_i = back_coord / (graph_seq_len + 1)
            back_j = back_coord % (graph_seq_len + 1)
            # print("coming from: ", back_coord, back_i, back_j)

            # int node_id, int node_letter_int, int node_pos, int read_letter_int, int read_pos
            # print(back_i, back_j, i, j)
            if back_j == j:  # insertion
                alignment.add_alignment(j_node[j - 1], -1, j_pos[j - 1], read_as_int[i - 1], i - 1)
            elif back_i == i:  # deletion
                alignment.add_alignment(j_node[j - 1], all_seq_as_int[j], j_pos[j - 1], -1, i - 1)
            else:  # match or miss
                alignment.add_alignment(j_node[j - 1], all_seq_as_int[j - 1], j_pos[j - 1], read_as_int[i - 1], i - 1)

            i = back_i
            j = back_j
            coord = back_coord
            # if j == current[1]:  # insertion
            #     candidate_mappings = (node.id, i - 1, "_", seq_letter, node_str_idx)
            # elif i == current[0]:  # deletion
            #     candidate_mappings = (node.id, i - 1, column_letter, "_", node_str_idx)
            # else:  # match-miss
            #     candidate_mappings = (node.id, i-1, column_letter, seq_letter, node_str_idx)

        # print(alignment.info)
        alignment.output_paf(graph, output_file=None, stdout=True)
    # need to free the memory
    free(dp_table)
