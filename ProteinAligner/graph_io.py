import os
import sys
import logging
from ProteinAligner.Node import Node
import pdb


def write_gfa(graph, gfa_path, colored=False):
    """
    output the graph as a gfa file

    :param graph: a graph object
    :param gfa_path: the output gfa file name/path
    :param colored: If I want to output the classes too
    """
    # maybe a dictionary of nodes and their classes
    # then a function that reads these colors and make an upset plot or something
    nodes = graph.nodes
    if os.path.exists(gfa_path):
        logging.warning("overwriting {} file".format(gfa_path))

    # Assuming graphs here are made from MSA and don't have overlaps
    overlap = "0M"
    f = open(gfa_path, "w+")

    for node in nodes.values():
        if not colored:
            line = str("\t".join(("S", str(node.id), node.seq)))
        else:
            colors = "|".join(list(node.colors))
            line = str("\t".join(("S", str(node.id), node.seq, colors)))
        f.write(line + "\n")

        for child in node.out_nodes:
            edge = str("\t".join(("L", str(node.id), "+", str(child),
                       "+", overlap)))
            f.write(edge + "\n")

    if len(graph.paths) != 0:
        for p_name in graph.paths.keys():
            segment = []
            for n in graph.paths[p_name]:
                if graph.nodes[n].id in graph.nodes:
                    segment.append(graph.nodes[n].id)
            segment = ",".join([str(x) + "+" for x in segment])
            overlaps = "0M," * (len(graph.paths[p_name]) - 1)
            overlaps = overlaps[0:len(overlaps) - 1]
            out_path = "\t".join(["P", p_name, segment, overlaps])
            f.write(out_path + "\n")

    f.close()


def read_gfa(gfa_path, colored=False):
    """
    read a gfa file and store the information in nodes dictionary

    :param gfa_path: path to the graph file to read
    :param k: the k value used for building k-mers
    :param colored: if True, then colors for each node are read and stored
    """
    if not os.path.exists(gfa_path):
        print("Error happened, aborting!!! Please check log file")
        logging.error("The file {} does not exist".format(gfa_path))
        sys.exit(1)

    nodes = dict()
    edges = []
    with open(gfa_path, "r") as in_graph:
        for line in in_graph:

            if line.startswith("S"):  # a node (Segment)
                line = line.strip().split()

                # 1 would be the identifier and 2 is the sequence
                ident = int(line[1])
                nodes[ident] = Node(ident, line[2])
                if colored:
                    nodes[ident].colors = set(line[3].split(":"))

            elif line[0] == "L":  # an edge
                edges.append(line)

        for e in edges:
            # assuming all forward nodes here because directed aa graph (no reverse complement)
            # and the graph represent the protein (from beginning to end of MSA)
            e = e.strip().split()
            parent = int(e[1])
            child = int(e[3])
            if (parent in nodes) and (child in nodes):
                nodes[parent].out_nodes.add(child)
                nodes[child].in_nodes.add(parent)

    return nodes
