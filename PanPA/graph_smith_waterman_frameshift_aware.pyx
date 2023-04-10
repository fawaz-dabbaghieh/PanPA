# distutils: language=c++

from PanPA.Graph cimport Graph
from PanPA.Alignment cimport Alignment
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free
from libcpp.string cimport string
from PanPA.constants import all_linear_sub_matrices
from PanPA.reverse_complement_fast import reverse_complement
from PanPA.constants import translation_table


cdef void print_dp_table(int all_seq_len, int read_len, int *dp_table, vector[int] graph_seq,
                         vector[int] read_seq):
    """
    This function prints the dp table for debugging purposes
    turns the 1D table back into 2D and prints each line
    """
    line = []
    cdef int i, row, j

    graph_seq_chr = []
    for i in graph_seq:
        graph_seq_chr.append(chr(i + 65))
    read_seq_chr = []
    for i in read_seq:
        read_seq_chr.append(chr(i + 65))

    print("      " + "  ".join(graph_seq_chr))
    for i in range(read_len + 1):
        row = i * (all_seq_len + 1)
        for j in range(all_seq_len + 1):
            line.append(dp_table[row + j])
        if i == 0:
            print(" ", line)
        else:
            print(read_seq_chr[i - 1], line)
        line = []
    print("\n")

cdef vector[string] align_to_graph_sw_fsa(Graph graph, str read, str read_name, bint print_dp,
                                      vector[int] sub_matrix, int gap_score, min_id_score) except *:
    """
    This function applies a Smith-Waterman algorithm but for graphs (Partial Order Alignment)
    The columns are constructed from the topologically ordered nodes
    The only difference is that instead of looking at j-1 for calculating the cell's score
    One needs to look at the last letter of all incoming nodes and to which j they belong to
    """
    # cdef vector[str] alignments
    cdef vector[string] alignments
    cdef string one_alignment
    alignments = []
    ########################################################

    """
    In this part, I am assigning pointers to structures in the graph object
    otherwise, if I access them with graph_obj.attribute, cython is adding extra python interaction lines
    that are checking whether the graph object has that attribute or not
    so just making a pointer avoids all these checks that might cause a slow down in the code in the main loop.
    """
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

    cdef vector[int] all_seq_as_int
    all_seq_as_int = graph.all_seqs
    #########################################################

    #########################################################
    """
    Defining some variables I need for the algorithm
    """
    cdef int which_node = -1  # counter to tell me in which node I am in the column loop

    # keeping track of all possible max scores in the whole dp_table
    cdef int global_max = 0
    cdef vector[int] global_max_coord  # cell coordinates of all max scores (if more than one exist)

    cdef int in_nodes_max, in_nodes_max_coord  # which in node gives the best score

    # self-explanatory variables
    cdef int n, dimensions, i, j, left_cell, current_cell, above_cell, diagonal_cell, row, l
    cdef int match_miss, deletion, insertion, local_max, previous_node, previous_j

    """
    for a certain column j, which node does that corresponds to and where in the node am I
    if the node has 4 letters, then if current_node_pos == 0 then I need to look at incoming nodes
    because I am at position 0 in the node sequence (left side of the node)
    otherwise I only need to look at j-1 when calculating the score
    """
    cdef int current_node_id, current_node_pos
    cdef int counter = 0  # for testing
    #########################################################

    #########################################################
    """
    Converting the read to integers
    declaring the dp table and filling it with 0 values

    and filling traceback table
    """
    cdef vector[int] read_as_int
    cdef int *dp_table  # dp table is an array with allocated memory
    cdef int read_len, graph_seq_len

    read_len = len(read)
    graph_seq_len = graph.all_seqs.size()

    dimensions = (read_len + 1) * (graph_seq_len + 1)
    dp_table = <int *> malloc(dimensions * sizeof(int))  # allocating memory

    # initializing the vector with 0s, otherwise random values will be filled
    for i in range(dimensions):
        dp_table[i] = 0

    """
    65 because for ASCII this is where the capital letters start so instead of A 65 it's now 0
    I do this because I linearized my substitution matrix, and need to use the ints to do fast access
    to this linearized substitution matrix when retrieving the pair-score
    """
    for c in read:
        read_as_int.push_back(ord(c) - 65)

    # a traceback table also filled with 0s, later will be filled with the index of where the score came from
    cdef int *traceback_table
    traceback_table = <int *> malloc(dimensions * sizeof(int))
    for i in range(dimensions):
        traceback_table[i] = 0

    #########################################################
    # main loop for the alignment

    for i in range(read_len):
        # adding one to keep first row and first column 0
        i += 1

        which_node = -1  # resetting the node counter

        # because the dp table is 1d, I need to convert between i,j to 1d coordinates
        row = i * (graph_seq_len + 1)
        for j in range(graph_seq_len):

            # counter += 1
            maximum = 0
            j += 1
            # print("i {}, j {}".format(i, j))
            """
            I have 3 situations:

            1- I am at the left border of a node and the node doesn't have in nodes (no more tracebacks)

            2- I am not at the left border of the node so I only need to look at i-1 j-1 in the dp-table

            3- I am at the left border and the node has in nodes
            This one then need the column jumping, because I need to look at the previous_j for each in node, which
            corresponds to the end of that node. e.g. node1:SVT-->node2:RPP and if I am R then I need to look at where
            T is in the dp-table because it's not necessarily at j-1 but j' that corresponds to T of node1 is definitely
            smaller than j because of the topological sorting
            """
            in_nodes.clear()
            in_nodes_size = 0

            current_node_id = j_node[j - 1]  # which node corresponds to the current j
            current_node_pos = j_pos[j - 1]  # where in that node am I
            if current_node_pos == 0:  # we are at a new node
                which_node += 1

            for l in range(intervals[which_node], intervals[which_node + 1]):
                in_nodes.push_back(in_nodes_vec[l])
            # in_nodes = in_nodes_map[which_node]
            in_nodes_size = in_nodes.size()

            # print(f"i {i}, j {j}, current_node {current_node_id}, current_node_pos {current_node_pos}, in_nodes {in_nodes}")

            # 1 and 2 situation
            if (current_node_pos != 0) or (current_node_pos == 0 and in_nodes_size == 0):
                # if (current_node_pos != 0) or :

                current_cell = row + j
                left_cell = current_cell - 1
                above_cell = current_cell - graph_seq_len - 1
                diagonal_cell = above_cell - 1

                """
                because the substitution matrix is linearized over all the alphabet + the character for the stop codon (27 characters)
                both the read and the sequences in the graph are integers with A being 0 and Z 25
                therefore, accessing a certain "cell" in the linearized 2d substitution matrix corresponds to
                first_letter * 27 + second_letter

                This adds a bit speed up compared to accessing key:value maps
                """
                match_miss = dp_table[diagonal_cell] + sub_matrix[read_as_int[i - 1] * 27 + all_seq_as_int[j - 1]]

                deletion = dp_table[left_cell] + gap_score
                insertion = dp_table[above_cell] + gap_score

                local_max = max(match_miss, deletion, insertion, 0)

                if local_max == match_miss:
                    traceback_table[current_cell] = diagonal_cell
                elif local_max == deletion:
                    traceback_table[current_cell] = left_cell
                elif local_max == insertion:
                    traceback_table[current_cell] = above_cell


                dp_table[current_cell] = local_max


            # situation 3
            elif (current_node_pos == 0) and (in_nodes_size != 0):

                # if in_nodes_size != 0:

                # for in nodes I only need to choose one so coord is an integer and not a vector like the global max
                in_nodes_max = 0
                in_nodes_max_coord = 0

                """
                minor optimizations
                current_cell, above_cell and insertion don't depend on the incoming edge but only on current j
                therefore, can only need to be calculated once before looping through all incoming edges
                """
                current_cell = row + j
                above_cell = current_cell - graph_seq_len - 1

                insertion = dp_table[above_cell] + gap_score

                for previous_node in in_nodes:

                    # node_ends tells me where to jump in the dp table
                    # so instead of j-1 I get some other previous_j
                    previous_j = node_ends[previous_node]

                    left_cell = row + previous_j
                    diagonal_cell = left_cell - graph_seq_len - 1

                    match_miss = dp_table[diagonal_cell] + sub_matrix[read_as_int[i - 1] * 27 + all_seq_as_int[j - 1]]
                    #
                    # print(i, j, read[i - 1], node.seq[current_node_pos],
                    #       sub_matrix[read_as_int[i-1] + all_seq_as_int[j-1]], in_nodes_size)

                    deletion = dp_table[left_cell] + gap_score

                    # updating in_nodes_max which will tell me which jump had the best score
                    local_max = max(match_miss, deletion, insertion, 0)
                    # print(f"2 read/graph {chr(read_as_int[i - 1] + 65)}/{chr(all_seq_as_int[j - 1] + 65)} match or miss score {match_miss}, deletion {deletion}, insertion {insertion}, local_max {local_max}")

                    if local_max > in_nodes_max:
                        # filling the traceback table with the potential max
                        if local_max == match_miss:
                            in_nodes_max_coord = diagonal_cell
                        elif local_max == deletion:
                            in_nodes_max_coord = left_cell
                        elif local_max == insertion:
                            in_nodes_max_coord = above_cell

                        in_nodes_max = local_max
                    # print(i, j, current_cell, left_cell, above_cell, diagonal_cell, local_max)
                # updating the global max and its coordinates
                # if in_nodes_max > global_max:
                #     global_max = local_max
                #     global_max_coord.clear()
                #     global_max_coord.push_back(current_cell)
                # elif in_nodes_max == global_max:
                #     global_max_coord.push_back(current_cell)

                traceback_table[current_cell] = in_nodes_max_coord
                dp_table[current_cell] = in_nodes_max

    # print(counter)
    if print_dp:
        print("This is the dp table\n")
        print_dp_table(graph_seq_len, read_len, dp_table, all_seq_as_int, read_as_int)
        print("This is the traceback table\n")
        print_dp_table(graph_seq_len, read_len, traceback_table, all_seq_as_int, read_as_int)

    ###############################################################################################
    # finding the max scores and coordinates in the DP table to traceback these best alignments
    for i in range(dimensions):
        if dp_table[i] > global_max:
            global_max = dp_table[i]
            global_max_coord.clear()
            global_max_coord.push_back(i)
        elif dp_table[i] == global_max:
            global_max_coord.push_back(i)

    # print(f"The max location in matrix is {global_max_coord} and the score is {global_max}")
    # traceback now
    cdef int coord, back_coord
    cdef Alignment alignment
    cdef int traceback_max, back_j, back_i, alignment_len
    cdef int test = 0
    cdef float alignment_score
    # print(j_node, j_pos)

    for coord in global_max_coord:
        alignment = Alignment(read_name, read_len, global_max)

        traceback_max = dp_table[coord]
        i = coord / (graph_seq_len + 1)
        j = coord % (graph_seq_len + 1)

        # while not first row or column, and not a 0 score (maybe alignment just ended in the middle)
        # while ((i != 0) and (j != 0)) and (dp_table[coord] != 0):
        while ((i != 0) and (j != 0)) and (dp_table[coord] != 0):
            # print("current: coord, i, j", coord, i, j)

            back_coord = traceback_table[coord]
            back_i = back_coord / (graph_seq_len + 1)
            back_j = back_coord % (graph_seq_len + 1)

            if back_j == j:  # insertion
                # print("insertion")
                alignment.add_alignment(j_node[j - 1], -1, j_pos[j - 1], read_as_int[i - 1], i - 1)
            elif back_i == i:  # deletion
                # print("deletion")
                alignment.add_alignment(j_node[j - 1], all_seq_as_int[j - 1], j_pos[j - 1], -1, i - 1)
            else:  # match or miss
                # print("match or mismatch")
                alignment.add_alignment(j_node[j - 1], all_seq_as_int[j - 1], j_pos[j - 1], read_as_int[i - 1], i - 1)

            i = back_i
            j = back_j
            coord = back_coord

        alignment.prepare_gaf(graph, min_id_score)
        alignment_score = alignment.id_score
        # print(f"the alignment score is {alignment_score} and the min id is {min_id_score}")
        if alignment_score >= min_id_score:
            alignments.push_back(alignment.gaf.encode())

    # need to free the memory, otherwise major memory leaks
    free(dp_table)
    free(traceback_table)
    return alignments


def align_dna_seq(gfa_file, seq, seq_name, sub_matrix_name, gap_s, min_id_score, print_dp):
    cdef Graph graph
    cdef vector[int] sub_matrix
    # for the codon translation I can generate a static vector[vector[int]] once and I just give it to this functions
    # which gives it to the alignment function to use
    # so I only construct it once when I run the program, i.e. construct it in _main and give it to this function later

    # for testing I can copy some of these files and build them in place with a simple setup script
    # like this python3 setup.py build_ext --inplace
    # or make a new subcommand for now and have the test there and keep building PanPA I guess
    for i in all_linear_sub_matrices["blosum62"]:
        sub_matrix.push_back(i)

    graph = Graph(gfa_file)
    alignments = align_to_graph_sw_fsa(graph, seq, seq_name, print_dp, sub_matrix, gap_s, min_id_score)