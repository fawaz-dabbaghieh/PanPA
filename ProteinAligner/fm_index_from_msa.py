import shellinford
import pickle
import os
import sys
import logging
from ProteinAligner.read_fasta import read_fasta

"""
The initial idea here is that this will be called with either all MSA's or the GFAs generated from the MSAs
Then for each graph, a concatenated string of all paths is constructed.

A list then for each graph is given to shellinford FM-Index which saves it into a binary file. When we have reads,
k-mers of minimizers are taken from the reads and ran through the fm.search function.

The order of the list of MSA' or GFAs need to be maintained (maybe a simple pickled dic in the same directory as the 
FM index) which tells us which doc_id corresponds to which graph or MSA.
"""


def build_from_msas(list_of_msas, index_location):
    """
    takes a list of MSA files and builds the FM-Index for all of them and a pickled dict that points
    to which file has which id

    :param list_of_msas: a list of fasta file paths
    :param index_location: where to save the FM index and the pickled dict
    """

    conc_sequences = dict()
    logging.info("Reading and concatenating sequences before building the FM-Index")
    for msa in list_of_msas:
        sequences = read_fasta(msa)
        conc_sequences[msa] = ""
        for s in sequences.values:
            conc_sequences[msa] += s.replace("-", "")  # removing gaps from MSA sequences

    msa_order = dict()
    idx = 0
    for msa in conc_sequences.keys():
        msa_order[idx] = msa
        idx += 1

    logging.info("Building the FM index")
    pickled_dict_loc = os.path.join(index_location, "msa_order.pickle")
    pickled_dict = open(pickled_dict_loc, "wb")
    pickle.dump(msa_order, pickled_dict)
    pickled_dict.close()
    fm_index_loc = os.path.join(index_location, "fm_index_from_msa.fm")
    shellinford.FMIndex(list(conc_sequences.values()), fm_index_loc)


def search_seed(seeds, index_location="."):
    """
    takes a seed and checked the fm index and returns to which MSA it should match to

    :param seeds: a list of seeds from one read or contig
    :param index_location: Directory where to keep the fm index
    """
    index_file = os.path.join(index_location, os.sep, "fm_index_all_msa.fm")
    pickled_loc = os.path.join(index_location, os.sep, "msa_order.pickle")

    if not os.path.exists(index_location):
        print("Error! Please check the log file")
        logging.error(f"The FM index file {index_file} was not found")
        sys.exit(1)

    if not os.path.exists(pickled_loc):
        print("Error! Please check the log file")
        logging.error(f"The order of MSAs pickled dict file {pickled_loc} was not found")
        sys.exit(1)

    fm = shellinford.FMIndex(filename=index_file)
    with open(pickled_loc, "rb") as handle:
        msa_order = pickle.load(handle)

    matches = dict()
    # for each seed I look at all the matches and count the matches for each msa (graph)
    for seed in seeds:
        for s in fm.search(seed):
            if msa_order[s.doc_id] in matches:
                matches[msa_order[s.doc_id]] += 1
            else:
                msa_order[s.doc_id] = 1

    return matches
