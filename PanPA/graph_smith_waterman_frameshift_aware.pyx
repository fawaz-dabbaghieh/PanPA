# distutils: language=c++
# cython: profile=True

from PanPA.Graph cimport Graph
from PanPA.Alignment cimport Alignment
from libcpp.vector cimport vector
from libc.stdlib cimport malloc, free
from libcpp.string cimport string

cdef void print_dp_table(int graph_seq_len, int read_len, int *dp_table, vector[int] graph_seq,
                         vector[int] read_seq):
    """
    This function prints the dp table for debugging purposes
    turns the 1D table back into 2D and prints each line
    """
    line = []
    cdef int i, row, j
    process_nucleotides = {0: "A", 1: "C", 2: "T", 3: "G", 4: "N"}

    out_dp = open("dp_table.csv", "w")

    graph_seq_chr = []
    for i in graph_seq:
        graph_seq_chr.append(chr(i + 65))

    read_seq_chr = []
    for i in read_seq:
        read_seq_chr.append(process_nucleotides[i])

    print("    " + " ".join(graph_seq_chr))
    out_dp.write("  " + " ".join(graph_seq_chr))
    out_dp.write("\n")
    for i in range(read_len + 2):
        row = i * (graph_seq_len + 1)
        for j in range(graph_seq_len + 1):
            line.append(dp_table[row + j])
        if i == 0 or i == 1:
            print(f"  {line}".replace("[", "").replace("]", "").replace(",", ""))
            out_dp.write(f" {line}".replace("[", "").replace("]", "").replace(",", ""))
            out_dp.write("\n")
        else:
            print(f"{read_seq_chr[i - 2]} {line}".replace("[", "").replace("]", "").replace(",", ""))
            out_dp.write(f"{read_seq_chr[i - 2]} {line}".replace("[", "").replace("]", "").replace(",", ""))
            out_dp.write("\n")
        line = []
    print("\n")


cdef vector[string] align_to_graph_sw_fsa(Graph graph, str read, str read_name, bint print_dp,
                                          vector[int] sub_matrix, int gap_score, int fs_score, float min_id_score,
                                          vector[vector[vector[int]]] codon_translate_0base,
                                          vector[int] neucleotide_to_int) except *:
    """
    This function applies a Smith-Waterman algorithm but for graphs (Partial Order Alignment)
    The columns are constructed from the topologically ordered nodes
    The only difference is that instead of looking at j-1 for calculating the cell's score
    One needs to look at the last letter of all incoming nodes and to which j they belong to
    """

    # cdef vector[str] alignments
    cdef vector[string] alignments
    cdef string one_alignment
    # alignments = []
    ########################################################
    cdef int codon
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
    # todo move this to outside and give as argument, so it's not initialized everytime the function is called

    cdef int nuc
    read_len = len(read)
    graph_seq_len = graph.all_seqs.size()

    # + 2 because of the 2 extra initialization lines for the jumps
    dimensions = (read_len + 2) * (graph_seq_len + 1)
    dp_table = <int *> malloc(dimensions * sizeof(int))  # allocating memory

    # initializing the vector with 0s, otherwise random values will be filled
    for i in range(dimensions):
        dp_table[i] = 0

    """
    65 because for ASCII this is where the capital letters start so instead of A 65 it's now 0
    I do this because I linearized my substitution matrix, and need to use the ints to do fast access
    to this linearized substitution matrix when retrieving the pair-score, much faster than using
    a hash table or dictionary for getting the scores
    """
    # instead of using ord(c) - 65 I need to use the simple map of ACTGN 01234
    # process_nucleotides = {"A":0, "C":1, "T":2, "G":3, "N":4}
    # for c in read:
    #     read_as_int.push_back(ord(c) - 65)
    for c in read:
        # print(f"converting letter {c} with ord {ord(c)}")
        nuc = ord(c)
        read_as_int.push_back(neucleotide_to_int[nuc])


    # a traceback table also filled with 0s, later will be filled with the index of where the score came from
    cdef int *traceback_table
    traceback_table = <int *> malloc(dimensions * sizeof(int))
    for i in range(dimensions):
        traceback_table[i] = 0

    #########################################################
    # main loop for the alignment
    # print(f"graph seq len is {graph_seq_len}")
    # -2 because I skip the first two nucleotide and consider the codon from its last nucleotide
    for i in range(read_len - 2):
        # adding one to keep first row and first column 0
        i += 2
        codon = codon_translate_0base[read_as_int[i - 2]][read_as_int[i - 1]][read_as_int[i]]
        i += 2
        which_node = -1  # resetting the node counter

        # because the dp table is 1d, I need to convert between i,j to 1d coordinates
        row = i * (graph_seq_len + 1)
        for j in range(graph_seq_len):

            # counter += 1
            # print(f"the value of j is {j}")
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
            in_nodes_size = in_nodes.size()

            # print(f"i {i}, j {j}, current_node {current_node_id}, current_node_pos {current_node_pos}, in_nodes {in_nodes}")
            # 1 and 2 situation
            if (current_node_pos != 0) or (current_node_pos == 0 and in_nodes_size == 0):
                # if (current_node_pos != 0) or :

                current_cell = row + j
                left_cell = current_cell - 1
                above_cell = current_cell - (graph_seq_len + 1) * 3
                diagonal_cell = above_cell - 1

                # framshift jumps
                diagonal_insertion = current_cell - (graph_seq_len + 1) * 4 - 1
                diagonal_deletion = current_cell - (graph_seq_len + 1) * 2 - 1

                """
                because the substitution matrix is linearized over all the alphabet + the character for the stop codon (27 characters)
                both the read and the sequences in the graph are integers with A being 0 and Z 25
                therefore, accessing a certain "cell" in the linearized 2d substitution matrix corresponds to
                first_letter * 27 + second_letter

                This adds a bit speed up compared to accessing key:value maps
                """
                # print(f"we are at row {i} and column {j}")
                # print(f"We are at cell {current_cell} with coords of {i},{j} with left_cell {left_cell}, above_cell {above_cell} and diagonal is {diagonal_cell}, fs deletion is {diagonal_deletion} and fs insertion is {diagonal_insertion}")
                # print(f"We are at codon {chr(codon + 65)} and target is {chr(all_seq_as_int[j-1] + 65)}, and the score is {sub_matrix[codon * 27 + all_seq_as_int[j - 1]]}")
                if codon == all_seq_as_int[j - 1]:
                    match_miss = dp_table[diagonal_cell] + 2
                else:
                    match_miss = dp_table[diagonal_cell] + 0
                # match_miss = dp_table[diagonal_cell] + sub_matrix[codon * 27 + all_seq_as_int[j - 1]]

                deletion = dp_table[left_cell] + gap_score
                insertion = dp_table[above_cell] + gap_score
                fs_deletion = dp_table[diagonal_deletion] + fs_score
                fs_insertion = dp_table[diagonal_insertion] + fs_score

                local_max = max(match_miss, deletion, insertion, fs_deletion, fs_insertion, 0)

                if local_max == match_miss:
                    traceback_table[current_cell] = diagonal_cell
                elif local_max == deletion:
                    traceback_table[current_cell] = left_cell
                elif local_max == insertion:
                    traceback_table[current_cell] = above_cell
                elif local_max == fs_deletion:
                    # print(f"I am in diagonal deletion at {current_cell} and taking from {diagonal_deletion}")
                    traceback_table[current_cell] = diagonal_deletion
                elif local_max == fs_insertion:
                    # print(f"I am in diagonal insertion at {current_cell} and taking from {diagonal_deletion}")
                    traceback_table[current_cell] = diagonal_insertion

                # print(f"the local max is {local_max}")
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
                above_cell = current_cell - (graph_seq_len + 1) * 3

                insertion = dp_table[above_cell] + gap_score

                for previous_node in in_nodes:

                    # node_ends tells me where to jump in the dp table
                    # so instead of j-1 I get some other previous_j
                    previous_j = node_ends[previous_node]

                    left_cell = row + previous_j
                    diagonal_cell = left_cell - (graph_seq_len + 1) * 3
                    # framshift jumps

                    diagonal_insertion = left_cell - (graph_seq_len + 1) * 4
                    diagonal_deletion = left_cell - (graph_seq_len + 1) * 2

                    # print(f"we are at row {i} and column {j}")
                    # print(f"We are at codon {chr(codon + 65)} and target is {chr(all_seq_as_int[j-1] + 65)}, and the score is {sub_matrix[codon * 27 + all_seq_as_int[j - 1]]}")
                    if codon == all_seq_as_int[j - 1]:
                        match_miss = dp_table[diagonal_cell] + 2
                    else:
                        match_miss = dp_table[diagonal_cell] + 0
                    # match_miss = dp_table[diagonal_cell] + sub_matrix[codon * 27 + all_seq_as_int[j - 1]]
                    #
                    # print(i, j, read[i - 1], node.seq[current_node_pos],
                    #       sub_matrix[read_as_int[i-1] + all_seq_as_int[j-1]], in_nodes_size)

                    deletion = dp_table[left_cell] + gap_score
                    fs_deletion = dp_table[diagonal_deletion] + fs_score
                    fs_insertion = dp_table[diagonal_insertion] + fs_score

                    # updating in_nodes_max which will tell me which jump had the best score
                    local_max = max(match_miss, deletion, insertion, fs_deletion, fs_insertion, 0)
                    if local_max > in_nodes_max:
                        if local_max == match_miss:
                            in_nodes_max_coord = diagonal_cell
                        elif local_max == deletion:
                            in_nodes_max_coord = left_cell
                        elif local_max == insertion:
                            in_nodes_max_coord = above_cell
                        elif local_max == fs_deletion:
                            # print(f"I am in diagonal deletion at row {current_cell // (graph_seq_len + 1)} and column {current_cell % (graph_seq_len + 1)} taking from row {diagonal_deletion // (graph_seq_len + 1)} and column {diagonal_deletion % (graph_seq_len + 1)} with current score {local_max} and previous score {dp_table[diagonal_deletion]}")
                            in_nodes_max_coord = diagonal_deletion
                        elif local_max == fs_insertion:
                            # print(f"I am in diagonal insertion at row {current_cell // (graph_seq_len + 1)} and column {current_cell % (graph_seq_len + 1)} taking from row {diagonal_deletion // (graph_seq_len + 1)} and column {diagonal_deletion % (graph_seq_len + 1)} with current score {local_max} and previous score {dp_table[diagonal_deletion]}")
                            in_nodes_max_coord = diagonal_insertion
                        # print(f"{current_cell}, {local_max}, {in_nodes_max_coord}")

                        # local_max = max(match_miss, deletion, insertion, 0)
                        # print(f"2 read/graph {chr(read_as_int[i - 1] + 65)}/{chr(all_seq_as_int[j - 1] + 65)} match or miss score {match_miss}, deletion {deletion}, insertion {insertion}, local_max {local_max}")

                        # if local_max > in_nodes_max:
                        #     # filling the traceback table with the potential max
                        #     if local_max == match_miss:
                        #         in_nodes_max_coord = diagonal_cell
                        #     elif local_max == deletion:
                        #         in_nodes_max_coord = left_cell
                        #     elif local_max == insertion:
                        #         in_nodes_max_coord = above_cell

                        in_nodes_max = local_max
                    # print(f"We are at {current_cell} and the previous j is {previous_j}, the local_max is {local_max} and the coordinates of it is {in_nodes_max_coord}, diagonal is at {diagonal_cell}, and fs deletion {diagonal_deletion} and fs insertion {diagonal_insertion}")
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
        # print("This is the traceback table\n")
        # print_dp_table(graph_seq_len, read_len, traceback_table, all_seq_as_int, read_as_int)


    ###############################################################################################
    # finding the max scores and coordinates in the DP table to traceback these best alignments
    for i in range(dimensions):
        # had to add the > 0 check because on rare occasions, if the read being aligned really have nothing matching
        # the algorithm for some reason was getting stuck
        if dp_table[i] > 0:
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
    # alignment = []
    cdef int traceback_max, back_j, back_i, alignment_len
    cdef int test = 0
    cdef float alignment_score
    # print(j_node, j_pos)

    for coord in global_max_coord:
        alignment = Alignment(read_name, read_len, global_max)

        traceback_max = dp_table[coord]
        i = coord // (graph_seq_len + 1)
        j = coord % (graph_seq_len + 1)

        # while not first row or column, and not a 0 score (maybe alignment just ended in the middle)
        # while ((i != 0) and (j != 0)) and (dp_table[coord] != 0):
        while ((i != 0) and (j != 0)) and (dp_table[coord] != 0):
            # print("current: coord, i, j", coord, i, j)

            back_coord = traceback_table[coord]
            back_i = back_coord // (graph_seq_len + 1)
            back_j = back_coord % (graph_seq_len + 1)

            # print(f"I am at {i + 2}, {j + 2} and came from {back_i + 2}, {back_j + 2}")
            # print(f"and the graph letter is {chr(all_seq_as_int[j-1] + 65)}")
            # print(f"and the codon is {read[i-4] + read[i-3] + read[i-2]}")
            if back_j == j:  # insertion
                # print("insertion")
                # add - to target
                # read pos is i - 4 because I want to position to represent the beginning of the codon
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": "_", "read_pos": i - 2,
                     "read_str": read[i - 2], "type": 0})
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": "_", "read_pos": i - 3,
                     "read_str": read[i - 3], "type": 0})
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": "_", "read_pos": i - 4,
                     "read_str": read[i - 4], "type": 0})

                alignment.path.append(j_node[j - 1])
                alignment.n_indels += 3
                # alignment.add_alignment(j_node[j - 1], -1, j_pos[j - 1], read_as_int[i - 1], i - 1)

            elif back_i == i:  # deletion
                # skipping a whole codon in the read, so ---
                # print("deletion")
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": chr(all_seq_as_int[j - 1] + 65),
                     "read_pos": i - 2, "read_str": "---", "type": 1})
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": chr(all_seq_as_int[j - 1] + 65),
                     "read_pos": i - 3, "read_str": "---", "type": 1})
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": chr(all_seq_as_int[j - 1] + 65),
                     "read_pos": i - 4, "read_str": "---", "type": 1})

                alignment.path.append(j_node[j - 1])
                alignment.n_indels += 3

                # alignment.add_alignment(j_node[j - 1], all_seq_as_int[j - 1], j_pos[j - 1], -1, i - 1)
            elif back_i == i - 3:  # match or mismatch
                # if I am correct, the current letter is the end of the codon and is i-2, to get the codon, it's i-2-2, i-2-1, and i-2-0
                codon = codon_translate_0base[read_as_int[i - 4]][read_as_int[i - 3]][read_as_int[i - 2]]

                graph_character = chr(all_seq_as_int[j - 1] + 65)
                if all_seq_as_int[j - 1] == codon:
                    cigar_type = 2
                    alignment.n_matches += 3
                else:
                    cigar_type = 3
                    alignment.n_mismatches += 3
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": graph_character, "read_pos": i - 2,
                     "read_str": read[i - 2], "type": cigar_type})
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": graph_character, "read_pos": i - 3,
                     "read_str": read[i - 3], "type": cigar_type})
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": graph_character, "read_pos": i - 4,
                     "read_str": read[i - 4], "type": cigar_type})

                alignment.path.append(j_node[j - 1])
                # add both target letter and
                # I need to take the codon at i-2, i-1, i
                # print("match or mismatch")
                # alignment.add_alignment(j_node[j - 1], all_seq_as_int[j - 1], j_pos[j - 1], read_as_int[i - 1], i - 1)

            elif back_i == i - 2:  # deletion frameshift
                # print(f"current {i}, {j} and backtrace are {back_i}, {back_j}")
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": chr(all_seq_as_int[j - 1] + 65),
                     "read_pos": i - 2, "read_str": read[i - 2], "type": 3})
                alignment.info.append(
                    {"node_id": j_node[j - 1], "node_pos": j_pos[j - 1], "node_str": chr(all_seq_as_int[j - 1] + 65),
                     "read_pos": i - 3, "read_str": read[i - 3], "type": 3})

                alignment.path.append(j_node[j - 1])
                alignment.n_mismatches += 2

            elif back_i == i - 4:  # insertion frameshift
                # todo I think I still need to check the insertion-match model, not just the match-insertion
                codon1 = codon_translate_0base[read_as_int[i - 5]][read_as_int[i - 4]][read_as_int[i - 3]]
                codon2 = codon_translate_0base[read_as_int[i - 4]][read_as_int[i - 3]][read_as_int[i - 2]]
                graph_character = chr(all_seq_as_int[j - 1] + 65)

                # hacky, but saves me time instead of re-writing a bunch of things in the Alignment class
                if all_seq_as_int[j - 1] == codon1:
                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 2,
                                           "read_str": read[i - 2],
                                           "type": 0})

                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 3,
                                           "read_str": read[i - 3],
                                           "type": 2})
                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 4,
                                           "read_str": read[i - 4],
                                           "type": 2})
                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 5,
                                           "read_str": read[i - 5],
                                           "type": 2})

                    alignment.n_matches += 3
                elif all_seq_as_int[j - 1] == codon2:

                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 2,
                                           "read_str": read[i - 2],
                                           "type": 2})
                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 3,
                                           "read_str": read[i - 3],
                                           "type": 2})
                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 4,
                                           "read_str": read[i - 4],
                                           "type": 2})

                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 5,
                                           "read_str": read[i - 5],
                                           "type": 0})
                    alignment.n_matches += 3
                else:
                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 2,
                                           "read_str": read[i - 2],
                                           "type": 0})

                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 3,
                                           "read_str": read[i - 3],
                                           "type": 3})
                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 4,
                                           "read_str": read[i - 4],
                                           "type": 3})
                    alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j - 1],
                                           "node_str": chr(all_seq_as_int[j - 1] + 65), "read_pos": i - 5,
                                           "read_str": read[i - 5],
                                           "type": 3})

                    alignment.n_mismatches += 3
                # alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j-1], "node_str": graph_character,"read_pos": i - 4, "read_str": read[i-4] + read[i-3] + read[i - 2], "cigar": cigar})

                # alignment.info.append({"node_id": j_node[j - 1], "node_pos": j_pos[j-1], "node_str": chr(all_seq_as_int[j-1] + 65), "read_pos": i - 4, "read_str": read[i-5] + read[i-4] + read[i-3] + read[i - 2], "cigar": cigar})
                alignment.path.append(j_node[j - 1])
                alignment.n_indels += 1

            i = back_i
            j = back_j
            coord = back_coord

        alignment.prepare_aa_gaf(graph)
        alignment_score = alignment.id_score
        if alignment_score >= min_id_score:
            alignments.push_back(alignment.gaf.encode())
    # need to free the memory, otherwise major memory leaks
    free(dp_table)
    free(traceback_table)
    return alignments
