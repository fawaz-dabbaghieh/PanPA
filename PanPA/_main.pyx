# distutils: language=c++
import sys
import os
import logging
import pickle
import multiprocessing as mp
from PanPA.reverse_complement_fast import reverse_complement
from libcpp.vector cimport vector
from libcpp.string cimport string
from PanPA.Graph cimport Graph
from PanPA.read_fasta import read_fasta, read_fasta_gen
from PanPA.graph_from_msa import msa_graph
from PanPA.index_sequences import *
from PanPA.graph_smith_waterman cimport align_to_graph_sw
from PanPA.graph_smith_waterman_frameshift_aware cimport align_to_graph_sw_fsa
from PanPA.constants import all_linear_sub_matrices


def exit_error(message):
    print("Something went wrong, please check the log file for more information")
    logging.error(message)
    sys.exit(1)


def check_output_file(output_file):
    if os.path.exists(output_file):
        exit_error(f"The file {output_file} already exists, not going to overwrite it")
    else:
        try:
            test = open(output_file, "w")
            test.close()
        except PermissionError:
            exit_error(f"Was not able to create {output_file}, permission denied")
        except MemoryError as e:
            exit_error(f"There was a memory error {e}")


def check_existance(files):
    """
    check if the paths of given files exist or not, if not, it will write a warning to the log file and remove
    these non-existent files from the list so the tool doesn't break because one or two wrong paths
    """
    not_exist = []
    for f in files:
        if not os.path.exists(f):
            not_exist.append(f)
            logging.warning(f"File {f} does not exist, will be skipped")

    for f in not_exist:
        files.remove(f)


def align_to_graph(seqs_dict, graphs, graph_files, sub_matrix, gap_score, fs_score,
                   seed_index, index_info, seed_limit, min_id_score, queue, printdp, dna):
    """
    Takes a dict with several reads and their names, and align them and add them to the queue
    """
    cdef Graph graph
    cdef bint print_dp
    cdef vector[string] alignments
    cdef int gap_s = gap_score
    cdef string a
    cdef int fs_s = fs_score
    # this maps a codon to an amino acid with ord(A)-65
    # ACTG are mapped with neucleotide_to_int to integers 01234
    # so for a codon like AGT we use codon_translate_0base[0][3][2] and we'll get an int
    # to turn the int back to an amino acid, e.g. chr(codon_translate_0base[0][3][2]+65)
    # and we'll get the amino acid S
    cdef vector[vector[vector[int]]] codon_translate_0base = [
        [[10, 13, 13, 10, 23], [19, 19, 19, 19, 23], [8, 8, 8, 12, 23], [17, 18, 18, 17, 23], [23, 23, 23, 23, 23]],
        [[16, 7, 7, 16, 23], [15, 15, 15, 15, 23], [11, 11, 11, 11, 23], [17, 17, 17, 17, 23], [23, 23, 23, 23, 23]],
        [[26, 24, 24, 26, 23], [18, 18, 18, 18, 23], [11, 5, 5, 11, 23], [26, 2, 2, 22, 23], [23, 23, 23, 23, 23]],
        [[4, 3, 3, 4, 23], [0, 0, 0, 0, 23], [21, 21, 21, 21, 23], [6, 6, 6, 6, 23], [23, 23, 23, 23, 23]],
        [[23, 23, 23, 23, 23], [23, 23, 23, 23, 23], [23, 23, 23, 23, 23], [23, 23, 23, 23, 23], [23, 23, 23, 23, 23]]]
    cdef vector[int] neucleotide_to_int = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
                                           4, 0, 0, 0, 0, 0, 2]

    # cdef vector[string] alignments_vec
    print_dp = printdp
    # if DNA is true, then instead of going through each seq I go through each seq's frames and add
    # that frame to the alignment output with RF:i:value (values are 1,2,3,-1,-2,-3)
    # frames_to_bytes = {1:b"1", 2:b"2", 3:b"3", -1:b"-1", -2:b"-2", -3:b"-3"}
    for s_name, s in seqs_dict.items():
        current_seq = dict()
        if dna:
            reverse = reverse_complement(s)
            current_seq = {0:[s_name, s], 1:[s_name, reverse]}
        else:
            current_seq = {0:[s_name, s]}
        #     reading_frames = translate(seq)
        # else:
        #     reading_frames = {0:seq}  # just a dict with one sequences

        # for frame, seq in reading_frames.items():
        for frame, content in current_seq.items():
            # print(f"I am at frame {frame} with content {content}")
            seq_name = content[0]
            seq = content[1]
            if dna:
                # todo need to think whether to exclude the forward or the reverse
                #   because in theory, one of them must be the correct one
                matches = query_dna_sequence(seq, seed_index, index_info)
                # print(matches)
            else:
                matches = query_aa_sequence(seq, seed_index, index_info)
            # print(matches)
            if 0 < seed_limit < len(matches):
                # keeping only the seeds up to the limit
                matches = matches[0:seed_limit]

            # print(f"The matches are {matches}")

            if matches:
                for i in matches:
                    if graphs[i] is None:
                        logging.warning("There's a seed without a matching graph")
                        continue

                    graph = graphs[i]
                    #print(f"aligning to graph {graph.name}")
                    # print(f"going to align {seq_name} to {graph.name}")
                    if dna:
                        alignments = align_to_graph_sw_fsa(graph, seq, seq_name, print_dp, sub_matrix,
                                                           gap_s, fs_s, min_id_score, codon_translate_0base,
                                                           neucleotide_to_int)
                    else:
                        alignments = align_to_graph_sw(graph, seq, seq_name, print_dp, sub_matrix,
                                                       gap_s, min_id_score)
                    # print(f"aligned {seq_name} to graph {i} and got {alignments.size()} alignments back")
                    # if I wanted to add another tag, I can do it here with the frame translation
                    # I mean add it for each alignment
                    # for a in alignments:
                    #     queue.put(a)
                    for a in alignments:
                        if not dna:  # it's a protein sequence
                            queue.put(a)
                        else:  # it's a DNA and could be reversed
                            if frame == 0:
                                queue.put(a)
                                # queue.put(a + b"\tst:Z:forward")
                            else:
                                # todo hacky and needs fixing later
                                a = a.replace(b"\t+\t", b"\t-\t")
                                queue.put(a)
                                # queue.put(a + b"\tst:Z:reverse")
                    # if not reading_frame:
                    #     queue.put(a)
                    # else:
                    #     queue.put(a + b"\tRF:i:" + reading_frame)
                    # I am not adding alignments to queue immediately to decrease the amount of synchronization needed
                    # while several processes are trying to add to the queue, after getting all alignments
                    # I add them to the queue
                    # for a in alignments:
                    #     alignments_vec.push_back(a)
                    # print(f"It took {time.perf_counter() - start} to align {seq_name} with length {len(seq)}")

                # for a in alignments_vec:
                #     queue.put(a)
    queue.put(b'0')  # sentinel


def align(graphs, index, graph_files, sub_matrix, args):
    seed_index = index['seed_index']
    index_info = index['index_info']
    files_index = index['files_index']

    counter = 0
    processes = []
    seq_counter = 0
    # tmp_out = 0
    queue = mp.Queue()
    out_gaf = open(args.out_gaf, "w")
    seqs_dicts = [dict()]
    dna = args.is_dna  # true or false
    for seq_name, seq in read_fasta_gen(args.in_seqs):
        # if DNA I can convert here and continue normally
        # the only difference is that I need a a reading frame tag
        # I can add the reading frame to the end of seq_name
        # and extract it inside align_to_graph function and add it to the alignment produced
        # so I don't really need to change much in the code

        # todo I need to add a counter for how many alignments done
        #   I can add a new queue maybe where each process returns an int of how many alignments it did
        #   I can sum them all and print the number
        #   but then it won't be consistent (not regular intervals)
        if seq_counter != 1000:  # preparing batch
            seqs_dicts[-1][seq_name] = seq
            seq_counter += 1

        else:
            seqs_dicts[-1][seq_name] = seq
            seq_counter = 0

            if len(seqs_dicts) != args.n_cores:  # we don't have enough batches yet
                seqs_dicts.append(dict())
            else:  # we have enough batches to deploy
                for seqs in seqs_dicts:  # for each batch
                    # (seqs_dict, graphs, graph_files, sub_matrix,
                    #  gap_score, seed_index, index_info, seed_limit, min_id_score, queue, printdp,
                    #  dna)
                    # I can also add the reading frame now, which is a binary string
                    # def align_to_graph(seqs_dict, graphs, graph_files, sub_matrix, gap_score, fs_score,
                    #                    seed_index, index_info, seed_limit, min_id_score, queue, printdp, dna):
                    p = mp.Process(target=align_to_graph, args=(seqs, graphs, graph_files,
                                                                sub_matrix, args.gap_score, args.fs_score, seed_index, index_info,
                                                                args.seed_limit, args.min_id_score, queue, False, dna,))
                    processes.append(p)

                for p in processes:
                    p.start()

                n_sentinals = 0
                while n_sentinals != args.n_cores:
                    a = queue.get()
                    if a == b'0':
                        n_sentinals += 1
                    else:
                        out_gaf.write(a.decode() + "\n")

                for p in processes:
                    p.join()

                processes = []
                queue = mp.Queue()
                n_sentinals = 0
                seq_counter = 0
                seqs_dicts = [dict()]

        # counter += 1
        # if counter % 200 == 0:
        #     logging.info(f"Processed {counter} reads so far...")

    # leftovers
    for seqs in seqs_dicts:
        processes.append(mp.Process(target=align_to_graph, args=(seqs, graphs, graph_files,
                                                                 sub_matrix, args.gap_score, args.fs_score, seed_index, index_info,
                                                                args.seed_limit, args.min_id_score, queue, False, dna,)))
    new_sent_len = len(processes)
    for p in processes:
        p.start()
    n_sentinals = 0
    while n_sentinals != new_sent_len:
        a = queue.get()
        if a == b'0':
            n_sentinals += 1
        else:
            out_gaf.write(a.decode() + "\n")

    for p in processes:
        p.join()

    out_gaf.close()
    logging.info("Done!")


def load_graph(graph_file, graph_n, graphs, lock):
    """
    Loads a graph file to a graph object and adds it to the queue
    """
    cdef Graph graph
    graph = Graph(gfa_file=graph_file)
    print(f"acquiring the lock for graph {graph_file}")
    lock.acquire()
    # out_graphs.append(graph)
    try:
        graphs[graph_n] = graph
    finally:
        print(f"releasing the lock for graph {graph_file}")
        lock.release()


def call_align_single(in_graph, seqs_dict, is_dna, sub_matrix, gap_score, fs_score, min_id_score, queue):
    cdef Graph graph
    cdef bint print_dp = False
    cdef vector[string] alignments
    cdef int gap_s = gap_score
    cdef string a
    cdef int fs_s = fs_score
    # this maps a codon to an amino acid with ord(A)-65
    # ACTG are mapped with neucleotide_to_int to integers 01234
    # so for a codon like AGT we use codon_translate_0base[0][3][2] and we'll get an int
    # to turn the int back to an amino acid, e.g. chr(codon_translate_0base[0][3][2]+65)
    # and we'll get the amino acid S
    cdef vector[vector[vector[int]]] codon_translate_0base = [
        [[10, 13, 13, 10, 23], [19, 19, 19, 19, 23], [8, 8, 8, 12, 23], [17, 18, 18, 17, 23], [23, 23, 23, 23, 23]],
        [[16, 7, 7, 16, 23], [15, 15, 15, 15, 23], [11, 11, 11, 11, 23], [17, 17, 17, 17, 23], [23, 23, 23, 23, 23]],
        [[26, 24, 24, 26, 23], [18, 18, 18, 18, 23], [11, 5, 5, 11, 23], [26, 2, 2, 22, 23], [23, 23, 23, 23, 23]],
        [[4, 3, 3, 4, 23], [0, 0, 0, 0, 23], [21, 21, 21, 21, 23], [6, 6, 6, 6, 23], [23, 23, 23, 23, 23]],
        [[23, 23, 23, 23, 23], [23, 23, 23, 23, 23], [23, 23, 23, 23, 23], [23, 23, 23, 23, 23], [23, 23, 23, 23, 23]]]
    cdef vector[int] neucleotide_to_int = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
                                           4, 0, 0, 0, 0, 0, 2]
    graph = in_graph

    for s_name, s in seqs_dict.items():
        current_seq = dict()
        if is_dna:
            reverse = reverse_complement(s)
            current_seq = {0:[s_name, s], 1:[s_name, reverse]}
        else:
            current_seq = {0:[s_name, s]}
        for frame, content in current_seq.items():
            # print(f"going to align {content}")
            seq_name = content[0]
            seq = content[1]
            if is_dna:
                # print(f"sending {seq} to graph aligner")
                alignments = align_to_graph_sw_fsa(graph, seq, seq_name, print_dp, sub_matrix,
                                                   gap_s, fs_s, min_id_score, codon_translate_0base,
                                                   neucleotide_to_int)
                # print(f"Alignments we got were {alignments}")
            else:
                alignments = align_to_graph_sw(graph, seq, seq_name, print_dp, sub_matrix, gap_score, min_id_score)

            for a in alignments:
                if not is_dna:  # it's a protein sequence
                    queue.put(a)
                else:  # it's a DNA and could be reversed
                    if frame == 0:
                        # print(f"trying to add {a} to the queue")
                        queue.put(a)
                        #queue.put(a + b"\tDNA:Z:forward")
                        # print(f"added {a} to the queue")
                    else:
                        # print(f"trying to add {a} to the queue")
                        a = a.replace(b"\t+\t", b"\t-\t")
                        queue.put(a)
                        # queue.put(a + b"\tRF:Z:reverse")
                        # print(f"added {a} to the queue")

    queue.put(b'0')


def align_single_graph(in_graph, in_seqs, is_dna, n_cores, sub_matrix, out_gaf,
                       gap_score, fs_score, min_id_score):
    """
    This function aligns directly the sequences in FASTA to a graph file given
    without the need to look at an index, but only aligns to one single target GFA
    """

    cdef Graph graph
    cdef vector[string] alignments
    batch_size = 1000
    graph = Graph(in_graph)


    processes = []
    seq_counter = 0
    queue = mp.Queue()
    out_gaf = open(out_gaf, "w")
    seqs_dicts = [dict()]
    for seq_name, seq in read_fasta_gen(in_seqs):
        if seq_counter != batch_size:  # preparing batch of sequences
            seqs_dicts[-1][seq_name] = seq
            seq_counter += 1

        else:
            seqs_dicts[-1][seq_name] = seq
            seq_counter = 0

            if len(seqs_dicts) != n_cores:  # we don't have enough batches yet
                seqs_dicts.append(dict())
            else:  # we have enough batches to deploy
                for seqs in seqs_dicts:  # for each batch
                    # def call_align_single(in_graph, seqs_dict, is_dna, sub_matrix, gap_score, fs_score, min_id_score, queue):

                    p = mp.Process(target=call_align_single, args=(graph, seqs, is_dna, sub_matrix, gap_score,
                                                                   fs_score, min_id_score, queue, ))
                    processes.append(p)

                for p in processes:
                    p.start()

                n_sentinals = 0
                while n_sentinals != n_cores:
                    a = queue.get()
                    if a == b'0':
                        n_sentinals += 1
                    else:
                        out_gaf.write(a.decode() + "\n")

                for p in processes:
                    p.join()

                processes = []
                queue = mp.Queue()
                n_sentinals = 0
                seq_counter = 0
                seqs_dicts = [dict()]

    for seqs in seqs_dicts:  # leftovers

        processes.append(mp.Process(target=call_align_single, args=(graph, seqs, is_dna, sub_matrix, gap_score,
                                                                    fs_score, min_id_score, queue,)))

    new_sent_len = len(processes)
    for p in processes:
        p.start()
    n_sentinals = 0
    while n_sentinals != new_sent_len:
        a = queue.get()
        if a == b'0':
            n_sentinals += 1
        else:
            out_gaf.write(a.decode() + "\n")

    for p in processes:
        p.join()

    out_gaf.close()
    logging.info("Done!")


def make_graph(in_msa, out_directory):
    """
    Generate a GFA file from an MSA and write it to the output directory given
    """
    cdef Graph graph
    sequences = read_fasta(in_msa)
    graph = msa_graph(sequences)

    # graph_name = ".".join(in_msa.split(os.sep)[-1].split(".")[:-1])
    graph_name = in_msa.split(os.sep)[-1] + ".gfa"
    graph_loc = os.path.join(out_directory, graph_name)
    graph.write_gfa(gfa_path=graph_loc)


def process_inputs(args, subcommand):
    """
    Processes the common input the 3 commands have, so checks if the input files are correct. Or if the list of files
    provided has correct file paths, or if the directory given as input is available and have files to process
    """

    # checking input directory
    input_files = []
    if args.in_files is None and args.in_list is None and args.in_dir is None:
        print("Error! Please check the log file")
        logging.error("You did not provide any in files, or a list of files, or an input directory.")
        sys.exit(1)

    # checking if input files are given as several arguments
    if args.in_files is not None:
        check_existance(args.in_files)
        if len(args.in_files) == 0:
            print("Error! Please check the log file")
            logging.error("There were no more input files to check after skipping the non-existent ones")
            sys.exit(1)

        input_files = args.in_files
    # checking if a list of file in a txt file was given and checking the files inside
    elif args.in_list is not None:
        if os.path.exists(args.in_list):
            with open(args.in_list) as in_file:
                for l in in_file:
                    input_files.append(l.strip())
        else:
            print("Error! Please check the log file")
            logging.error(f"The file {args.in_list} given does not exist")
            sys.exit(1)

        check_existance(input_files)
        if len(input_files) == 0:
            print("Error! Please check the log file")
            logging.error("There were no more input files to check after skipping the non-existent ones")
            sys.exit(1)

    # checking if a directory was given
    elif args.in_dir:
        if not os.path.exists(args.in_dir):
            print("Error! Please check the log file")
            logging.error(f"The input directory {args.in_dir} you gave does not exist")
            sys.exit(1)

        if subcommand in {"build_index", "build_gfa"}:
            endings = {"fasta", "fa", "gz", "fas", "fna"}
        else:  # then the third subcommand is to align to graphs
            endings = {"gfa"}

        for f in os.listdir(args.in_dir):
            if f.split(".")[-1] in endings:
                input_files.append(os.path.join(args.in_dir, f))
            else:
                logging.warning(f"the file {f} has an unrecognized extension, allowed are [fasta, fa, gz, fas, fna]")

    input_files.sort()
    return input_files


def _main(sys_argv, args, msa_name=None):
    cdef Graph graph
    cdef vector[int] sub_matrix

    if args.subcommands is None:
        print("Please provide a subcommand, check -h --help for help.")
        sys.exit(1)

    ########################################################## Building index (indexing) ##############################
    if args.subcommands == "build_index":

        if args.seeding_alg not in {"wk_min", "k_mers"}:
            print("Error! Please check the log file")
            logging.error(f"The seeding algorithm {args.seeding_alg} not recognized")
            sys.exit(1)

        # checking if the size of the window makes sense
        if args.seeding_alg == "wk_min":
            if args.window_size == 1:
                args.seeding_alg = "k_mers"
                logging.warning("Window size of 1 results in using k_mers indexing")

        if os.path.exists(args.out_index):
            print("Error! Please check the log file")
            logging.error(f"The file {args.out_index} already exists, give another path or name")
            sys.exit()

        input_msas = process_inputs(args, args.subcommands)

        # if 0 given, then the limit is the number of MSAs which is the upper limit
        # if a positive integer given, then it's limited to number of MSAs in case it was bigger
        # negative values are not accepted
        if args.seed_limit > 0:
            # a seed can be at most present in all graphs
            if args.seed_limit > len(input_msas):
                seed_membership_limit = len(input_msas)
            else:
                seed_membership_limit = int(args.seed_limit)
        elif args.seed_limit == 0:
            seed_membership_limit = len(input_msas)
        else:
            print("Error! Please check the log file")
            logging.error("The seed limit has to be 0 or a positive integer")
            sys.exit(1)

        seed_index = dict()
        index_info = dict()
        files_index = dict()
        # creating a graph from each msa and putting them in a dict
        for idx, msa in enumerate(input_msas):
            logging.info(f"Reading the sequences of {msa}")
            sequences = read_fasta(msa)

            logging.info(f"Indexing {len(sequences)} sequences from {msa} and using {args.seeding_alg} seeding method")

            tmp = msa.split(os.sep)[-1] + ".gfa"
            # tmp = ".".join(msa.split(os.sep)[-1].split(".")[:-1])
            files_index[tmp] = idx  # dict of file names without extension and pos
            # indexing these sequences
            index_sequences(seed_index, sequences, idx, k=args.k_mer_size,
                            w=args.window_size, seeding_alg=args.seeding_alg)
            logging.info(f"Finished indexing {msa} and size of index is {len(seed_index)}")

        logging.info(f"Sorting the index and limiting it to {seed_membership_limit} matches per seed")
        sort_index_and_limit(seed_index, seed_membership_limit)
            # logging.info(f"Generating a DAG from the {msa}")

        # maybe hacky, but hiding the seeding type, k_mer length and window length in the seed_index dict
        if args.seeding_alg == "k_mers":
            logging.info(f"The seeding method is k-mers with a size of {args.k_mer_size}")
            index_info["type"] = "k_mers"
            index_info["ksize"] = args.k_mer_size
            index_info["wsize"] = 0

        else:
            logging.info(f"The seeding method is wk minimizers with k-mer size of {args.k_mer_size} and "
                         f"window size of {args.window_size}")
            index_info["type"] = "wk_min"
            index_info["ksize"] = args.k_mer_size
            index_info["wsize"] = args.window_size

        logging.info(f"Writing the index as a binary file to disk to {args.out_index}")
        try:
            out_index = open(args.out_index, "wb")
        except Exception as e:
            print("Error! Please check the log file")
            logging.error(f"When opening f{args.out_index} this happene: {e}")
            sys.exit(1)

        # the index is these 3 dictionaries
        to_pickle = {"seed_index":seed_index, "index_info":index_info, "files_index":files_index}
        pickle.dump(to_pickle, out_index)
        out_index.close()
        logging.info("Done!")

    ########################################################## Building GFA #####################################
    if args.subcommands == "build_gfa":

        input_msas = process_inputs(args, args.subcommands)
        # checking if the output directory already exists or not
        # and weather the user have permission to actually make a new directory
        # or there's space to make a new one
        if not os.path.exists(args.out_dir):
            try:
                logging.info(f"The output directory {args.out_dir} does not exist, trying to create one")
                os.mkdir(args.out_dir)
            except PermissionError:
                print("Error! Please check the log file")
                logging.error(f"Was not able to create the directory {args.out_dir} due to permission error")
                sys.exit(1)

        # generating output directory
        # gfa_output = os.path.join(args.out_dir, "graphs")
        # try:
        #     os.mkdir(gfa_output)
        # except FileExistsError:
        #     print("Error! Please check the log file")
        #     logging.error(f"The directory {gfa_output} already exists, change it or give another output directory")
        #     sys.exit(1)
        # except PermissionError:
        #     print("Error! Please check the log file")
        #     logging.error(f"A PermissionError occurred when making the directory {gfa_output}, check your permissions")
        #     sys.exit(1)

        # upper limit to number of cores that can be used
        if args.n_cores > os.cpu_count():
            args.n_cores = os.cpu_count()

        logging.info(f"Beginning to reads in MSAs and build graphs to {args.out_dir}")
        counter = 0
        processes = []
        checkpoint = int(len(input_msas)/10)
        if checkpoint == 0:  # avoid divide by 0
            checkpoint = 1

        for idx, msa in enumerate(input_msas):
            process = mp.Process(target=make_graph, args=(msa, args.out_dir,))
            processes.append(process)
            counter += 1
            if len(processes) == args.n_cores:
                for p in processes:
                    p.start()
                for p in processes:
                    p.join()
                # emptying to prepare the next batch of graphs
                processes = []

            if counter % checkpoint == 0:
                logging.info(f"So far {counter} MSAs have been processed")

        # processing the leftovers
        if processes:
            for p in processes:
                p.start()
            for p in processes:
                p.join()

        logging.info("Done!")
    ########################################################## Aligning #####################################
    if args.subcommands == "align":

        # printing available substitution matrices for the user
        if args.sub_matrix_list:
            available_matrices = []
            for k in all_linear_sub_matrices.keys():
                available_matrices.append(k)
            message = f"The following substitution matrices are available: {available_matrices}"
            print(message)
            sys.exit(0)

        check_output_file(args.out_gaf)

        if not 0 <= args.min_id_score <= 1:
            exit_error("The min alignment identity score needs to be between 0 and 1")

        # checking the substitution matrix the user chose exists or not
        if not args.sub_matrix in all_linear_sub_matrices:
            available_matrices = []
            for k in all_linear_sub_matrices.keys():
                available_matrices.append(k)
            message = f"The substitution matrix {args.sub_matrix} is not present in the constants file\n" \
                      f"The substitution matrix {args.sub_matrix} is not present in the constants file\n" \
                      f"The following substitution matrices are available: {available_matrices}"
            exit_error(message)
        else:
            for i in all_linear_sub_matrices["blosum62"]:
                sub_matrix.push_back(i)

        # checking and loading the index
        if args.index is None:
            exit_error("You did not provide an index file, if you don't have one, please run the build_index subcommand"
                       " to build an index from the MSAs")

        if not os.path.exists(args.index):
            exit_error(f"The file {args.index} provided does not exist")
            sys.exit(1)

        else:
            index_file = open(args.index, "rb")
            index = pickle.load(index_file)

        if args.seed_limit < 0:
            exit_error(f"The seed limit you chose {args.seed_limit} cannot be smaller than 0")

        # process_inputs will return the graph files
        graph_files = process_inputs(args, args.subcommands)

        if args.n_cores > os.cpu_count():
            args.n_cores = os.cpu_count()
        graphs = [None] * len(index['files_index'])

        # generating a graph object for each gfa file and adding it to a list
        # sorted the files first to match the original indices
        logging.info("Reading the graphs and loading them into memory...")
        counter = 0
        # loading graphs in a multiprocess
        # graph_batch = []
        # checkpoint = int(len(graph_files)/10)
        # if checkpoint == 0:  # avoiding 0 division later
        #     checkpoint = 1
        # initialization
        """
        This used to have a queue and I add the loaded graph to a queue
        then I loop through the processes, extract the graph objects from the queue
        and add it to the list "graphs"
        but it kept hanging at .join() and I couldn't find out why, so I decided to
        it might be that the graph objects are too big and the OS buffer is getting full
        or so I read online and it's hanging
        
        so I'll simply initialize a dictionary with keys being indices of range(n_graphs) 
        and I fill it by using lock and release on the threads
        
        did not work as well so decided to give up on multiprocessing loading
        need to check this again
        """

        # graphs = []
        for f in graph_files:
            logging.info(f"Loading graph {f}")
            graph = Graph(gfa_file=f)
            f_name = f.split(os.sep)[-1]
            # f_name = ".".join(f.split(os.sep)[-1].split(".")[:-1])
            idx = index["files_index"][f_name]
            graphs[idx] = graph
            # graphs.append(graph)
            logging.info(f"The graph {f} has {len(graph.nodes)} nodes and is loaded")


        logging.info(f"Aligning sequences to graphs and outputting to {args.out_gaf}")

        # I don't need a whole other function for align_dna
        # I can give a DNA flag, and the flag tells align_to_graph_ss
        # that this is a DNA sequence and I add the reading frame info the the seq name
        # so align_to_graph_ss extracts the reading frame name and adds it to the alignment
        # output if there was an alignment output
        align(graphs, index, graph_files, sub_matrix, args)
        # if not args.is_dna:
        #     # print("going into align_aa")
        #     align_aa(graphs, index, graph_files, sub_matrix, args)
        # else:
        #     align_dna(graphs, index, graph_files, sub_matrix, args)


    ########################################################## align single #####################################

    if args.subcommands == "align_single":
        # printing available substitution matrices for the user
        if args.sub_matrix_list:
            available_matrices = []
            for k in all_linear_sub_matrices.keys():
                available_matrices.append(k)
            message = f"The following substitution matrices are available: {available_matrices}"
            print(message)
            sys.exit(0)

        # repetitive code with align subcommand
        # But for now it's fine I guess
        if not 0 <= args.min_id_score <= 1:
            exit_error("The min alignment identity score needs to be between 0 and 1")

        # checking the substitution matrix the user chose exists or not
        if not args.sub_matrix in all_linear_sub_matrices:
            available_matrices = []
            for k in all_linear_sub_matrices.keys():
                available_matrices.append(k)
            message = f"The substitution matrix {args.sub_matrix} is not present in the constants file\n" \
                      f"The substitution matrix {args.sub_matrix} is not present in the constants file\n" \
                      f"The following substitution matrices are available: {available_matrices}"
            exit_error(message)
        else:
            for i in all_linear_sub_matrices["blosum62"]:
                sub_matrix.push_back(i)

        if args.n_cores > os.cpu_count():
            args.n_cores = os.cpu_count()

        if args.in_graph is None:
            exit_error("You did not provide a graph to align against")
        if not os.path.exists(args.in_graph):
            exit_error(f"The graph file {args.in_graph} does not exist")

        if args.in_seqs is None:
            exit_error("You did not provide a FASTA file with sequence to align against")
        if not os.path.exists(args.in_seqs):
            exit_error(f"The sequence file {args.in_seqs} does not exist")
        # checks if I can write to the output or not
        check_output_file(args.out_gaf)

        # takes all the args
        align_single_graph(args.in_graph, args.in_seqs, args.is_dna, args.n_cores, sub_matrix, args.out_gaf,
                           args.gap_score, args.fs_score, args.min_id_score)


    # TODO: I should add a stdout option to output alignments to stdout
    #     add a cutoff for alignment length as well, not just alignment identity
    #     i need to stop loading all graphs and only load the ones that have a hit in the index
    #     not as easy as I thought, as the loading will happen in a separate process so the same graph might get
    #     loaded several times which is not optimal

    """
    seems like the good way to do this is to start a p = multiprocessing.Process(target=alignment_fun, args(*args,)
    then I add this process to some list of processes, then I loop through all these processes in the list and do
    join on them, so this makes sure that I don't start the big loop (the one with the reads) before all the others
    finished already

    this might help me if I decide to write the Gotoh algorithm
    https://github.com/joachimwolff/algorithmsInBioinformatics/blob/master/Wolff_Presentations/Gotoh.pdf
    """
