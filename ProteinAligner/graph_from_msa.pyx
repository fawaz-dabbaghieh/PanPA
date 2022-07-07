# distutils: language=c++

import logging
from ProteinAligner.Node cimport Node
from ProteinAligner.Graph cimport Graph

"""
The idea here is simple
Let's take the following small example

01234567
-MEPTPEQ
---T--MA
MSETQSTQ

At column 0 we have two lists with the size := number of sequences
One is the previous columns (called previous_nodes) and current columns (current_nodes)
At first they are both empty, but then current is filled with [0,0,1]
where 1 is the key to the node constructed with sequence M

In the next iteration I sync previous and current. There's nothing to do here

Next, previous is [0,0,1] and current is [2,0,3] (2 is M and 3 is S from column 1)
when syncing, M will get connected to S. And current becomes previous

Let's say we are in column 3, we use the knowledge in already_seen to not make two T nodes
But only one with the information that T comes from both seq2 and seq3

In column 4, current is [8,0,9] and previous is [6,7,7] (T has id 7)
6 which is p gets connected to 8 which is T
T which 7 is connected to Q which is 9
We see that the middle T has a gap next, so we move it to current and current
become [8,7,9], the idea is that because of the gap, I keep pushing the previous letter or node
Until there's no more gap in column 6 where we have M now, and T is connected to M

So it's just previous and current columns following each other, and I keep connecting nodes
Or pushing the nodes where there's a gap until there isn't a gap.
"""


def sync_lists(nodes, current, previous):
    """
    syncing the previous column and current column lists

    The idea is that we have 4 possibilities for each row
    And 0 represent when we have a gap "-"
    1- previous has no node (has 0) and current is some amino acid
    nothing to do here

    2- previous position is a node and current is 0
    push the previous node to current replacing 0

    3- previous 0 and current 0
    nothing to do

    4- previous position is a node and current is a node
    connect them with an edge

    :param nodes: a dictionary of node ids and objects
    :param current: a list with nodes in the current column
    :param previous: a list with nodes in the previous column
    """
    cdef int i
    cdef Node current_node, previous_node
    # todo the same node shouldn't have a self loop (this is an old todo and not sure what to do about it)

    already_connected = set()

    for i in range(len(current)): # they both have the same size
        if (current[i] == 0) and (previous[i] != 0):
            # in case there was an amino acid in some sequence
            # then there's a gap, I push this amino acid (which is a node now)
            # to the next list, and I keep pushing it until the gap ends
            # so I can connect it with the next amino acid (node)
            current[i] = previous[i]

        # adding edges
        elif (current[i] != 0) and (previous[i] != 0):
            # node should be just a pointer to nodes[previous[i]]
            # so any changes to node should change nodes[previous[i]]
            # node = nodes[previous[i]]
            # node.add_child(nodes[current[i]])
            if not (current[i], previous[i]) in already_connected:  # todo this did not work
                previous_node = nodes[previous[i]]
                current_node = nodes[current[i]]

                previous_node.add_child(current_node)
                already_connected.add((current[i], previous[i]))

        else:
            pass
        # previous is current now and we have a new current
        # previous = current
        # current = [0] * len(previous)


def msa_graph(sequences):
    """
    build graph from MSA

    :param sequences: sequences dictionary
    :param seq_names: a dictionary of name of sequence in fasta file and potential shorter name
    """
    cdef Node node

    nodes = dict()
    previous_nodes = [0] * len(sequences)
    current_nodes = [0] * len(sequences)
    node_id = 1
    n_seqs = len(sequences[list(sequences.keys())[0]])
    seq_names_keys = list(sequences.keys())

    for j in range(n_seqs):

        # sync previous and current
        sync_lists(nodes, current_nodes, previous_nodes)
        previous_nodes = current_nodes
        current_nodes = [0] * len(previous_nodes)
        column = [x[j] for x in sequences.values()]
        already_seen = dict()

        for i, aa in enumerate(column):

            if aa != "-":  # skipping gaps
                if aa in already_seen.keys():
                    # adding the same amino acid at the same position
                    # this returns the last node created in this column with
                    # the same amino acid
                    node = nodes[already_seen[aa]]
                    color = seq_names_keys[i]
                    node.colors.add(color)
                    current_nodes[i] = node.identifier

                else:
                    # make a new node
                    node = Node(node_id, aa)
                    color = seq_names_keys[i]
                    node.colors.add(color)
                    nodes[node_id] = node
                    current_nodes[i] = node_id
                    # this tells me where I already saw this aa
                    # in which sequence (basically in which row)
                    already_seen[aa] = node_id
                    node_id += 1

    # syncing last column
    sync_lists(nodes, current_nodes, previous_nodes)

    graph = Graph(nodes=nodes)

    graph.compact()
    # graph.sort()
    # graph.add_paths()

    return graph
