# distutils: language=c++

import os
import logging
import sys
import json
import pickle
import collections
from PanPA.Node cimport Node


cdef class Graph:

    def __init__(self, gfa_file=None, paths=False, nodes=None):
        """
        If the gfa_file path is given, then the graph is initialized with that gfa graph
        Otherwise, the graph object will have emtpy containers

        To make the alignments fast, only vectors were used to store important information

        Using the following example graph, this will be explained:
        node id:   1   2  4   5    6   8   12  14
        node seq:  MS  M  E  PTPE  T  QST  MA  Q
        and the edges are 1->4, 2->4, 4->5, 4->6, 5->8, 5->14, 6->8, 6->12

        topological sort will be: MS  M  E  PTPE  T  QST  MA  Q

        Note that topological sort is not unique, exchanging MS and M at the beginning is still a valid sort


        After sorting, we concatenate all sequences to get the "text" we are aligning to (MSMEPTPETQSTMAQ), stored here
        as an integer vector with A->Z 0->25 in graph.all_seqs

        the vectors that are needed for the smith-waterman of this graph are as follows:

        j_pos is of len(all_seqs) and tells me for each j (column in DP-table) at which position am i in the node. e.g:
        the 0 1 2 3 corresponds to PTPE
               (M  S)(M)(E)(P  T  P  E) ....so on
        j_pos: [0, 1, 0, 0, 0, 1, 2, 3, 0, 0, 1, 2, 0, 1, 0]

        j_node is of len(all_seqs) tells for each j (column in DP-table) at which node id am I in. e.g.:
        5 5 5 5 corresponds to PTPE
        j_node: [1, 1, 2, 4, 5, 5, 5, 5, 6, 8, 8, 8, 12, 12, 14]

        so with j_pos and j_node I know for each j in which node I am and where in the node am I.


        topologically sorted nodes
        sorted: [1, 2, 4, 5, 6, 8, 12, 14]

        node_ends is same length as sorted (number of nodes), and tells me at what j the node ends. e.g: the number
        8 in node_ends corresponds to the column 8 in the DP-table and that's where the E is of PTPE
        node_ends: [2, 3, 4, 8, 9, 12, 14, 15]

        intervals is of length len(sorted) + 1, and it tells me the interval of in_nodes I need to look at, to get
        the incoming nodes of a certain node. e.g.: the first and second intervals are 0s because node 1 and 2 have no
        incoming edges, then we have the third interval 0-2 which is 0 and 1, looking at in_nodes[0] and in_nodes[1]
        gives back 0 and 1, taking sorted[0] and sorted[1] we get that node 4 has node 1 and node 2 as incoming edges

        now let's look a the final interval 6-8 which is 6 and 7, in_nodes[6] = 3 and in_nodes[7] = 5, then taking those
        and looking at sorted[3] and sorted[5] we get nodes 5 and 8, which are the incoming edges for node 14
        intervals: [0,0,0,2,3,4,5,6,8]
        in_nodes: [0,1,2,2,4,4,3,5]

        In the alignment algorithm everytime we have a 0 in j_pos, this means we entered a new node, and we can keep
        track at which node we are at, and how to look for its incoming edges.
        If we are at the left border of the node (at j_pos 0) then when calculating the score we need to look at all
        incoming nodes (intervals, in_nodes, sorted) tells us exactly who are these nodes and what previous j we need
        to jump to, if j_pos is not 0, then we only need to look at j-1
        """

        self.nodes = dict()
        self.paths = dict()
        # self.node_ends = dict()
        # self.column_pos = dict()
        # self.node_end = []
        # self.j_node = []
        # self.j_pos = []
        self.sorted = []
        self.seq_len = 0

        if gfa_file is not None:
            if os.path.exists(gfa_file):
                self.read_gfa(gfa_file, paths)
                self.name = gfa_file.split(os.sep)[-1]
                self.prepare_graph()

            else:
                print("Error! Check log file")
                logging.error("The file {} does not exist".format(gfa_file))
                sys.exit(1)

        # this will fill all the necessary vectors for aligning
        elif nodes is not None:
            self.nodes = nodes
            self.prepare_graph()
            self.name = ""
        else:  # nothing was giving to initialize the graph object
            self.nodes = dict()
            self.name = ""

    def __len__(self):
        return len(self.nodes)

    def __sizeof__(self):
        total_size = 0
        for n in self.nodes.values():
            total_size += n.__sizeof__()

        total_size += sys.getsizeof(self.sorted) + sys.getsizeof(self.all_seqs)  + sys.getsizeof(self.paths)
        total_size += sys.getsizeof(self.j_node)  + sys.getsizeof(self.j_pos)

        return total_size/1000000  # in mb

    cdef void prepare_graph(self) except *:
        """
        This function prepares the graph for alignment, calls the sorting function
        and fill the necessary containers that are needed by graph_align
        """
        # if the graph wasn't compacted originally, I don't want to compact it here
        # because the node ids will change then during compacting here then I will align to a compacted graph that will
        # slightly differ from the GFA given and this might cause a lot of confusion.
        # self.compact()
        if not self.sorted:
            self.sort()
        cdef Node node
        cdef int intervals_counter = 0
        cdef int i, n, j, seq_size

        cdef int idx = 0
        cdef int counter = 0

        nodes_to_sorted = dict()
        for idx, n in enumerate(self.sorted):
            nodes_to_sorted[n] = idx

        self.intervals.push_back(0)  # first value zero because I need to always take i and i-1 for the interval
        for n in self.sorted:
            node = self.nodes[n]
            seq_size = len(node.seq)
            counter += seq_size
            self.node_ends.push_back(counter)

            for i in range(seq_size):
                self.j_node.push_back(node.identifier)
                self.j_pos.push_back(i)

            for c in node.seq:
                self.all_seqs.push_back(ord(c) - 65)
            # for the one with intervals
            intervals_counter += node.in_nodes.size()
            self.intervals.push_back(intervals_counter)
            idx += 1

            for j in node.in_nodes:
                # self.sorted.index was very stupid and slow, I should've noticed this before
                # self.in_nodes.push_back(self.sorted.index(j))
                self.in_nodes.push_back(nodes_to_sorted[j])


    cdef void sort(self) except *:
        """
        Checks if graph is toplogically sorted, if not
        calls the toplogical sorting function
        """
        cdef int n

        if not self.sorted:
            self.top_sorting()

            if len(self.sorted) != len(self):  # It did not sort the graph
                print("Error! Check logs")
                logging.error("The graph was not sorted for some reason...")
                sys.exit(1)

            if not self.sorted:
                logging.error("The graph cannot be topologically sorted")
                sys.exit()
            elif len(self.sorted) != len(self.nodes):
                logging.error("The sorted list of nodes does not equal the number of nodes \n"
                              "Something went wrong, investigate please!")
                sys.exit()

    cdef void add_paths(self) except *:
        """
        Adds the different paths in the graph to the paths dictionary in the graph object
        """
        cdef int n
        cdef str color
        cdef Node node

        if not self.sorted:
            self.sort()

        if self.paths:
            self.paths.clear()

        for n in self.sorted:
            node = self.nodes[n]
            for color in node.colors:
                if color not in self.paths:
                    self.paths[color] = [n]
                else:
                    self.paths[color].append(n)

    cdef void nodes_info(self, output_file) except *:
        """
        writes out nodes information as JSON file

        :param output_file: path to output file
        """
        cdef dict nodes_info = dict()

        out_file = open(output_file, "w")
        for n in self.nodes.values():
            nodes_info[n.identifier] = {"id": n.identifier, "colors": list(n.colors), "seq": n.seq}
        out_file.write(json.dumps(nodes_info))
        out_file.close()

    cdef bint merge_end(self, int n) except *:
        """
        merges the end of a node if possible
        :param n: the node to try to merge its end
        """
        cdef Node node = self.nodes[n]
        cdef Node child = self.nodes[node.out_nodes[0]]
        cdef int cc_id = 0
        cdef Node childs_child
        cdef size_t one = 1
        cdef list out_nodes
        # if child has one parent, it can be merged otherwise not

        if (child.in_nodes.size() == one) and (node.colors == child.colors):
            '''
            adding the new children to n because
            n --> child --> child's_child1
                    \--> child's_child2
                    
            after merging n and child
            n+child --> child's child1
                   \--> child's child2
            '''
            node.out_nodes.clear()
            out_nodes = []
            for cc_id in child.out_nodes:
                out_nodes.append(cc_id)

            for cc_id in out_nodes:
                childs_child = self.nodes[cc_id]
                childs_child.remove_parent(child)
                node.add_child(childs_child)
                # self.nodes[n].out_nodes.add(new_out)

            # updating sequence of n assuming no overlaps to remove
            node.seq += child.seq

            # updating the information in children of child
            # in_nodes of child's children need to be updated
            # for nn in self.nodes[child].out_nodes:
            #     try:
            #         self.nodes[nn].in_nodes.remove(child)
            #     except KeyError:
            #         pdb.set_trace()
            #     self.nodes[nn].in_nodes.add(n)

            # remove the merged node
            del self.nodes[child.identifier]
            return True
        return False

    cdef void compact(self) except *:
        """
        compacts the graph
        """

        # I need to loop through a separate list because self.nodes will be changing
        cdef list list_of_nodes = list(self.nodes.keys())
        cdef int n
        cdef Node node
        for n in list_of_nodes:
            if n in self.nodes:
                while True:
                    # only one child (maybe can be merged)
                    node = self.nodes[n]
                    if node.out_nodes.size() == 1:
                        if not self.merge_end(n):
                            break
                    else:
                        break

        self.sorted.clear()
        self.add_paths()

    cdef void top_sorting(self) except *:
        """
        Performs Kahn's algorithm for topological sorting
        
        Kahn's algorithm needs to remove edges
        The easiest way was to make a deep copy of the graph
        If later I find this too problematic, I'll try the DFS based algorithm

        new_graph = copy.deepcopy(graph)
        For some reason I cannot call deepcopy on my graph
        there seems to be a problem with deep copy in python3
        it runs an infinite recursion and runs out of stack
        many stackoverflow posts talked about this problem

        so I'll just make a simple dictionary with the connections
        and use that to sort (less memory anyway than copying the whole graph)
        pseudo code from https://en.wikipedia.org/wiki/Topological_sorting
        
        """

        cdef dict nodes = dict()
        cdef Node n, new_node, child
        cdef int k, node_id
        cdef list sorted_nodes = []
        cdef list children = []
        cdef set starting_nodes = set()
        cdef bint empty

        # when using cython dictionaries, using items is faster than accessing only the key
        # for example then using that to get the value
        for k, n in self.nodes.items():
            new_node = Node(n.identifier)
            new_node.in_nodes = set([x for x in n.in_nodes])
            new_node.out_nodes = set([x for x in n.out_nodes])
            nodes[new_node.identifier] = new_node

        # sorted_nodes = []
        # starting_nodes = set()

        # getting starting nodes (nodes with no incoming edges)
        for n in nodes.values():
            empty = n.in_nodes.empty()
            if empty:
                starting_nodes.add(n.identifier)

        # Kahn's algorithm
        while starting_nodes:
            node_id = starting_nodes.pop()
            sorted_nodes.append(node_id)
            new_node = nodes[node_id]
            children = list(new_node.out_nodes)
            for child_id in children:
                # removing edge
                child = nodes[child_id]
                new_node.remove_child(child)
                # nodes[node_id].out_nodes.remove(child)
                # nodes[child].in_nodes.remove(node_id)

                # if child has no more parents
                # add to starting nodes set
                empty = child.in_nodes.empty()
                if empty:
                    starting_nodes.add(child.identifier)

        # if there is one or more edges left
        # then the graph is not DAG (or some other problem happened)
        for k, n in nodes.items():
            empty = n.in_nodes.empty()
            if not empty:
                self.sorted.clear()
        else:
            self.sorted = sorted_nodes


    cdef void read_gfa(self, str gfa_path, paths=True) except *:
        """
        read a gfa file and store the information in nodes dictionary

        :param gfa_path: path to the graph file to read
        :param paths: if True, then the path lines are read and added to the graph
        """

        cdef Node parent, child
        cdef list edges = []

        if not os.path.exists(gfa_path):
            print("Error happened, aborting!!! Please check log file")
            logging.error("The file {} does not exist".format(gfa_path))
            sys.exit(1)

        with open(gfa_path, "r") as in_graph:
            for line in in_graph:

                if line.startswith("S"):  # a node (Segment)
                    line = line.strip().split()

                    # 1 would be the identifier and 2 is the sequence
                    ident = int(line[1])
                    self.nodes[ident] = Node(ident, line[2])

                elif (line[0] == "P") and paths:
                    # path lines look like this
                    # P	path_14	11+,12-,13+	4M,5M
                    # I am ignoring the overlaps because I am assuming that protein graphs don't have overlaps
                    line = line.strip().split()
                    self.paths[line[1]] = [x[:-1] for x in line[2].split(",")]  # -1 to remove the +,-

                elif line[0] == "L":  # an edge
                    edges.append(line)

            for e in edges:
                # assuming all forward nodes here because directed aa graph (no reverse complement)
                # and the graph represent the protein (from beginning to end of MSA)
                e = e.strip().split()
                parent = self.nodes[int(e[1])]
                child = self.nodes[int(e[3])]
                if (parent.identifier in self.nodes) and (child.identifier in self.nodes):
                    # only add child because it a DAG
                    # and add_child will take care of filling in and out nodes
                    parent.add_child(child)
                    # self.nodes[parent].add_child(self.nodes[child])
                    # self.nodes[parent].out_nodes.add(child)
                    # self.nodes[child].in_nodes.add(parent)

    cdef void write_gfa(self, str gfa_path) except *:
        """
        output the graph as a gfa file

        :param gfa_path: the output gfa file name/path
        """

        cdef int k
        cdef Node node

        # maybe a dictionary of nodes and their classes
        # then a function that reads these colors and make an upset plot or something
        if os.path.exists(gfa_path):
            logging.warning("overwriting {} file".format(gfa_path))

        # Assuming graphs here are made from MSA and don't have overlaps
        overlap = "0M"
        f = open(gfa_path, "w+")
        # todo add an option to output the edge count when writing the GFA
        #    an edge count between two nodes is just how many shared colors between then
        for k, node in self.nodes.items():
            line = str("\t".join(("S", str(node.identifier), node.seq, f"LN:i:{len(node.seq)}" )))
            f.write(line + "\n")

            for child in node.out_nodes:
                edge = str("\t".join(("L", str(node.identifier), "+", str(child),
                           "+", overlap)))
                f.write(edge + "\n")

        if len(self.paths) != 0:
            for p_name, path in self.paths.items():
                # segment = []
                # for n in path:
                #     if n in self.nodes:
                #         segment.append(n)
                #
                # print(segment)
                segment = "+,".join([str(x) for x in path]) + "+"

                overlaps = ",".join(["0M"]*(len(path) - 1))
                out_path = "\t".join(["P", p_name, segment, overlaps])
                f.write(out_path + "\n")

        f.close()


    def pickle_info(self):
        cdef Node n
        in_nodes = dict()
        for n in self.nodes.values():
            in_nodes[n.identifier] = list(n.in_nodes)

        node_ends = self.node_ends

        j_pos = list(self.j_pos)
        j_node = list(self.j_node)
        all_seqs = list(self.all_seqs)

        all_info = dict()
        all_info["in_nodes"] = in_nodes
        all_info["node_ends"] = node_ends
        all_info["j_pos"] = j_pos
        all_info["j_node"] = j_node
        all_info["all_seqs"] = all_seqs

        with open('graph_info.pickle', 'wb') as handle:
            pickle.dump(all_info, handle)

    # cpdef void graph_align(self, str read,str read_name, bint print_dp, vector[int] sub_matrix, int gap_score,
    #                        output_file, min_id_score) except *:
    #
    #     if print_dp:
    #         print("j_pos:", self.j_pos)
    #         print("j_node:", self.j_node)
    #         print("sorted:", self.sorted)
    #         print("node_ends:", self.node_ends)
    #         print("in_nodes_to_sorted_index:", self.in_nodes_map)
    #         print("intervals:", self.intervals)
    #         print("in_nodes_conc:", self.in_nodes)
    #         print("\n\n")
    #     align_to_graph_sw(self, read, read_name, print_dp,
    #                       sub_matrix, gap_score, output_file, min_id_score)
    #
    # cpdef void graph_align_ed(self, str read, bint print_dp) except *:
    #     if print_dp:
    #         print("j_pos:", self.j_pos)
    #         print("j_node:", self.j_node)
    #         print("sorted:", self.sorted)
    #         print("node_ends:", self.node_ends)
    #         print("in_nodes_to_sorted_index:", self.in_nodes_map)
    #         print("intervals:", self.intervals)
    #         print("in_nodes_conc:", self.in_nodes)
    #         print("\n\n")
    #     align_to_graph_ed(self, read, print_dp)

    cdef str path_seq(self, list path):
        """
        Given a list of node ids, return a string of the path sequence
        or empty string if the path does not exist (there is no ordered walk through all these nodes)
        """
        cdef Node current_node, previous_node
        cdef list ordered_path = []
        cdef int i

        nodes_to_sorted = dict()
        for idx, n in enumerate(self.sorted):
            nodes_to_sorted[n] = idx

        NP = collections.namedtuple('NP', 'name position')
        for n in path:
            pos = nodes_to_sorted[n]
            ordered_path.append(NP(name=n, position=pos))

        # this way I sort the nodes in the path according to their position in the topological sort
        ordered_path.sort(key = lambda x: x.position)

        # now I check that each node is a child of the previous one, if so, then the path exists
        for i in range(1, len(ordered_path)):
            current_node = self.nodes[ordered_path[i].name]
            previous_node = self.nodes[ordered_path[i-1].name]
            if previous_node.is_parent_of(current_node):
                continue
            else:
                print("The ordered path is {}\n but node {} is not child of node {}".format(ordered_path,
                                                                                            current_node.identifier,
                                                                                            previous_node.identifier))
                print("The path is not a valid path")
                sys.exit()

        # AA(name=aa, count=count_list[aa])
        path_seq = ""
        for n in ordered_path:
            current_node = self.nodes[n.name]
            path_seq += current_node.seq

        return path_seq
