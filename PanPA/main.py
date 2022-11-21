#!/usr/bin/env python3

import sys
import argparse
import logging
from PanPA._main import _main

"""
1: Build Graph and Index
I want to build the FM-index and graphs from MSA input. I think the best way is to take a directory with FASTA
files inside, each FASTA is an MSA of some protein.
Then I build a GFA for each protein and also generate the FM-index.
However, there are a lot of similarities here, maybe it's a waste of time/space to concatenate all strings together
and generate the FM-index for all of them. Anyway, for now I'll do it this way and measure how much time it's adding
to each alignment.

2: Align
Takes the directory where the FM-index, the pickled dict that tells me which which doc corresponds to which MSA
these doc ids are related to the FM-index library I'm using, and takes the directory where all the GFA files are
Takes a FASTA file of reads to align, option for DNA or AA, if DNA, need to translate first into 6 different reading
frames, then take k-mer or minimizers from each read and run it by FM-index and see to which graph(s) there are hits
and align to these graphs.
"""

"""
Another way to do this, is maybe to just do it on the fly. So take the directory with all the MSAs, turn them into
graphs, also while reading the sequences I build a map of k-mers or minimizers and to which graph they belong to
and do the alignment.

Need to think about if a k-mer or a minimizer belongs to more than one graph, I guess I need to take a majority vote
so for one contig, check all the hits and see to which graphs they belong and the graph with the most hits I align to 
"""

def main():
    parser = argparse.ArgumentParser(description='Protein Graphs Aligner', add_help=True)
    subparsers = parser.add_subparsers(help='Available subcommands', dest="subcommands")

    parser._positionals.title = 'Subcommands'
    parser._optionals.title = 'Global Arguments'


    parser.add_argument("--log_file", dest="log_file", type=str, default="log.log",
                        help="The name/path of the log file. Default: log.log")

    parser.add_argument("--log_level", dest="log_level", type=str, default="INFO",
                        help="The logging level [DEBUG, INFO, WARNING, ERROR, CRITICAL]. Default: INFO")

    ############################################## indexing the MSAs ################################
    indexing_msas = subparsers.add_parser('build_index', help='Building an index from MSAs')

    indexing_msas.add_argument("-f", "--fasta_files", nargs="+", dest="in_files", type=str, default=None,
                                 help="Input MSA(s) in fasta format, one or more file space-separated")

    indexing_msas.add_argument("-l", "--fasta_list", dest="in_list", type=str, default=None,
                                 help="a text file with all input MSAs paths each on one new line")

    indexing_msas.add_argument("-d", "--in_dir", metavar="IN_DIR", dest="in_dir", default=None,
                                 type=str, help="Directory path containing one or more amino acid MSA in FASTA "
                                                "format (gzipped allowed)")

    indexing_msas.add_argument("-o", "--out_index", dest="out_index", default="index.pickle",
                                 type=str, help="The output index file name")

    indexing_msas.add_argument("--seeding_alg", dest="seeding_alg", default="k_mers",
                                 type=str, help="Seeding algorithm. Choices: k_mers, wk_min. Default: k_mers")

    indexing_msas.add_argument("-k", "--kmer_size", metavar="K-MER", dest="k_mer_size", default=5,
                                 type=int, help="K-mer size for indexing the sequencing. Default: 5")

    indexing_msas.add_argument("-w", "--window", metavar="WINDOW", dest="window_size", default=8,
                                 type=int, help="Window size when using w,k-minimizers instead of k-mers "
                                                "for indexing. Default:8")

    indexing_msas.add_argument("--seed_limit", dest="seed_limit", default=5,
                                 type=int, help="Indicates how many graphs can a seed belong to. Default: 5, "
                                                "give 0 for no limit")

    ############################################## building graph and indexes ################################
    building_graphs = subparsers.add_parser('build_gfa', help="building GFAs out of MSAs")


    building_graphs.add_argument("-f", "--fasta_files", nargs="+", dest="in_files", type=str, default=None,
                                 help="Input MSA(s) in fasta format, one or more file space-separated")

    building_graphs.add_argument("-l", "--fasta_list", dest="in_list", type=str, default=None,
                                 help="a text file with all input MSAs paths each on one new line")

    building_graphs.add_argument("-d", "--in_dir", metavar="IN_DIR", dest="in_dir", default=None,
                                 type=str, help="Directory path containing one or more amino acid MSA in FASTA "
                                                "format (gzipped allowed)")

    building_graphs.add_argument("-c", "--cores", metavar="CORES", dest="n_cores", default=1,
                          type=int, help="Numbers of cores to use for aligning")

    building_graphs.add_argument("-o", "--out_dir", metavar="OUT_DIR", dest="out_dir", default=".",
                                 type=str, help="Output directory where the index files and graphs from the MSAs are stored")

    ############################################## Aligning ################################

    aligning = subparsers.add_parser("align", help="aligning sequences given to graphs")

    aligning.add_argument("-g", "--gfa_files", nargs="+", dest="in_files", type=str, default=None,
                          help="Input GFA graphs, one or more file space-separated")

    aligning.add_argument("-l", "--gfa_list", dest="in_list", type=str, default=None,
                          help="a text file with all input graphs paths each on one new line")

    aligning.add_argument("-d", "--in_dir", metavar="GRAPHS", dest="in_dir", default=None,
                          type=str, help="Path to directory with GFA files")

    aligning.add_argument("--index", dest="index", default=None,
                          type=str, help="Path to pickled index file generated in the build step")

    aligning.add_argument("-r", "--seqs", metavar="SEQS", dest="in_seqs", default=None,
                          type=str, help="The input sequences to align in fasta format")

    aligning.add_argument("--dna", dest="is_dna", default=False, action="store_true",
                          help="Give this flag if the query sequences are DNA and not AA")

    aligning.add_argument("-c", "--cores", metavar="CORES", dest="n_cores", default=1,
                          type=int, help="Numbers of cores to use for aligning")

    aligning.add_argument("--sub_matrix", dest="sub_matrix", default="blosum62",
                                 type=str, help="Substitution matrix to use for alignment, default: blosum62")

    aligning.add_argument("--sub_matrix_list", dest="sub_matrix_list", default=False, action="store_true",
                          help="When given, a list of possible substitution matrices will be given")

    aligning.add_argument("-o", "--out_gaf", metavar="GAF", dest="out_gaf", default="alignments.gaf",
                          type=str, help="Output alignments file path")

    aligning.add_argument("--gap_score", dest="gap_score", default=-3, type=int,
                          help="The gap score to use for the alignment, default: -3")

    aligning.add_argument("--min_id_score", dest="min_id_score", default=0.7, type=float,
                          help="minimum alignment identity score for the alignment to be outputted, [0,1]")

    aligning.add_argument("--seed_limit", dest="seed_limit", default=3, type=int,
                          help="How many graphs can each seed from the query sequence have hits to, default: 3")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        print("You did not provide any arguments\n"
              "Try to use -h or --help for help\n")
        sys.exit()

    log_file = args.log_file

    logging.basicConfig(filename=log_file, filemode='w',
                        format='[%(asctime)s] %(message)s',
                        level=getattr(logging, args.log_level.upper()))

    logging.info(" ".join(["argument given:"] + sys.argv))

    _main(sys.argv, args)


if __name__ == "__main__":
    main()
