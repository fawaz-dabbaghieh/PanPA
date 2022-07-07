from collections import namedtuple, deque
import pickle
import sys
# from libcpp.deque cimport deque


def return_smallest(seq, k):
    """
    return smallest k-mer in the seq

    :param seq: the sequence
    :param k: the length of the kmer
    """
    cdef int len_seq
    cdef int k_size = k
    len_seq = len(seq)
    kmers = []

    for i in range(len_seq - k_size + 1):
        kmers.append(seq[i:i+k_size])

    return min(kmers)


def wk_minimizer(seq, k, w):
    """
    wk minimizer implementation

    :param seq: the string to find the minimizers for
    :param w: the window size
    :param k: the k-mer size
    """
    cdef int w_size = w
    cdef int k_size = k
    cdef int wk2 = w_size + k_size - 2
    cdef int wk = w_size + k_size
    cdef int len_seq = len(seq)
    cdef int i, j

    # minimizers = dict()

    # window smaller than sequence, no minimizers returned
    if (wk - 1) > len_seq:
        yield ""

    previous = ""
    for i in range(len_seq - wk2):
        sub_seq = seq[i: i + wk - 1]
        kmer = sub_seq[0:k_size]  # first kmer in window
        for j in range((w_size + k_size - 1) - k_size + 1):
            tmp = sub_seq[j:j+k_size]
            if tmp < kmer:
                kmer = tmp
        # min_kmer = return_smallest(seq[i: i + wk -1], k)
        # yield min_kmer
        if kmer != previous:
            previous = kmer
            yield kmer


def wk_minimizer_deque(seq, k, w):
    """
    This algorithm is taken from "Weighted minimizer sampling improves long read mapping" paper Algorithm 1
    this is a linear time wk-minimizer enumeration algorithm

    This keeps adding one more seed to each sequence at the beginning, the problem is in
    (if not minimizers) because if it stays empty, nothing will be added
    but there's no other way that simply adding the first seed
    I guess I should remove it later after making sure it is not actually the minimizer of the first window
    """

    cdef int len_seq = len(seq)
    cdef int k_size = k
    cdef int w_size = w
    cdef int i

    # cdef the_range = len_seq - k_size + 1
    minimizers = []
    my_deque = deque([])
    if len_seq < k_size + w_size - 1:
        return minimizers

    for i in range(len_seq - k_size + 1):

        current = (i, seq[i:i+k])
        # if the current seed is smaller than the last one in the deque, we pop the last one
        # we keep popping until current is not smaller anymore and we add it then to the end of the deque
        while my_deque and my_deque[-1][1] > current[1]:
            my_deque.pop()
        my_deque.append(current)
        # if the seeds to the left are outside of the range we pop them
        if my_deque[0][0] <= i - w:
            while my_deque[0][0] <= i - w:
                my_deque.popleft()  # discarding out of range k-mer

        # we don't have any minimizers yet we add one to the list
        if not minimizers:
            minimizers.append(my_deque[0][1])
        else:
            if minimizers[-1] != my_deque[0][1]:
                minimizers.append(my_deque[0][1])

    # This check is needed because I always add the first seed to get things started
    # in the if no minimizers. Otherwise minimizers will always stay empty and there's nothing to compare against
    # But now the first seed might be the minimal one, so I compare it to the next minimal one, if it's
    # bigger or equal, I can remove it, if it's smaller then it should stay
    if len(minimizers) > 1:
        if minimizers[0] <= minimizers[1]:
            del minimizers[0]
    else:
        return minimizers
    # hacky fix because in the first window all minimizers are included
    # return minimizers[w_size:] + [min(minimizers[0:w_size])]
    return minimizers


def extract_kmers(seq, k, *args):
    """
    A generator function that returns k-mers of the sequence given
    """
    cdef int i
    cdef int len_seq = len(seq)
    cdef int k_size = k

    if k_size > len_seq:
        yield ""
    elif k_size == len_seq:
        yield seq
    else:
        for i in range(len_seq - k_size + 1):
            yield seq[i:i+k_size]


def index_sequences(seed_index, sequences, graph_id, k, w, seeding_alg):
    """
    Takes a dictionary of sequences returned from reading the MSA
    And extracts seeds (k-mers or w-k-minimizers)
    """
    # probably hacky, but I didn't want to write the same code twice for k-mers and wk-minimizers
    # also didn't want to put the if k_mers inside the for seq loop because why should I ask this question everytime
    # but my functions are generators, so I need to have them in a loop
    # so now both functions take the same arguments, and I can put them in the seeding_funct variable, so I only
    # need to ask the questions whether to use k-mers or wk-minimizers once, and then I can use seeding_funct()

    if seeding_alg == "k_mers":
        seeding_funct = extract_kmers
    else:
        seeding_funct = wk_minimizer_deque

    # to normalize against number of sequences, instead of adding 1 every time I see this seed
    # I add the normalized count
    seed_count = 1/float(len(sequences))
    current_msa_seeds = dict()
    for seq in sequences.values():
        no_gaps_seq = seq.replace("-", "")  # removing - because it's an MSA
        for seed in seeding_funct(no_gaps_seq, k, w):
            # if the seed showed up 5 times and we have 4 sequences, then it will be 1
            # the count can get to more than 1 if there was a repeat for that seed
            if seed not in current_msa_seeds:
                current_msa_seeds[seed] = seed_count
            else:
                current_msa_seeds[seed] += seed_count

    for seed, counts in current_msa_seeds.items():
        if seed in seed_index:
            seed_index[seed].append([graph_id, counts])
        else:
            seed_index[seed] = [[graph_id, counts]]


def sort_index_and_limit(seed_index, seed_membership_limit):
    for key, value in seed_index.items():
        value.sort(key = lambda x : x[1], reverse=True)
        # updating each seed and only keeping the limit of how many graphs
        # ordered by their counts most to least
        seed_index[key] = [x[0] for x in value[0:seed_membership_limit]]



def query_sequence(seq, seed_index, index_info):
    """
    :param seq: a sequence that will be aligned
    :param index: the seed index dictionary
    """

    method = index_info["type"]
    seeds = []
    # extracting seeds
    if method == "wk_min":
        kmer = index_info["ksize"]
        window = index_info["wsize"]
        for seed in wk_minimizer_deque(seq, kmer, window):
            seeds.append(seed)
    else:
        kmer = index_info["ksize"]
        for seed in extract_kmers(seq, kmer):
            seeds.append(seed)

    # I count which graphs are mostly represented and order them based on the counts
    matches = dict()
    for seed in seeds:
        if seed in seed_index:
            for g in seed_index[seed]:
                if g not in matches:
                    matches[g] = 1
                else:
                    matches[g] += 1

    return sorted(matches, key=lambda k: matches[k], reverse=True)  # returns a sorted list
