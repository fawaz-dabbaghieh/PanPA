"""
This is an old version which is slow and bad, a linear version of wk minmizers is written in
index_sequences.pyx
"""
def return_smallest(seq, k):
    """
    return smallest k-mer in the seq

    :param seq: the sequence
    :param k: the length of the kmer
    """
    kmers = []

    for i in range(len(seq) - k + 1):
        kmers.append(seq[i:i+k])

    return min(kmers)


def wk_minimizer(string, w, k):
    """
    wk minimizer implementation

    :param string: the string to find the minimizers for
    :param w: the window size
    :param k: the k-mer size
    """

    minimizers = dict()
    if (w-1+k) > len(string):
        return minimizers

    for i in range(len(string) - (w + k - 2)):
        min_kmer = return_smallest(string[i:i+w+k-1], k)
        if min_kmer not in minimizers:
            minimizers[min_kmer] = 1
        else:
            minimizers[min_kmer] += 1

    return minimizers
