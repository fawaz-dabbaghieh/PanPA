# distutils: language=c++

import logging
from PanPA.constants import translation_table
from PanPA.reverse_complement_fast import reverse_complement

"""
Some notes from the Biopython library
NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid 
or a stop codon.  These are translated as "X".  Any invalid codon 
(e.g. "TA?" or "T-A") will throw an error.
"""
def inner_loop(seq):
    translations = [[], [], []]
    cdef int n, i, j
    n = len(seq)
    for i in range(0, n - n % 3, 3):

        for j in range(3):
            # the codon interval in the sequence
            start = i + j
            end = i + j + 3
            if end > n:  # end of sequence
                continue
            codon = seq[start:end].upper()

            # todo I can remove these conditions if I use my linearized translation table
            if "N" in codon:
                translations[j].append("X")  # we don't know what amino acid
            elif translation_table[codon] == "_":  # stop codon
                translations[j].append('[')
            elif codon not in translation_table:
                logging.warning("Warning! Codon {} not found in the translation "
                                "table and X is used for this codon".format(codon))
                translations[j].append("X")
            else:
                translations[j].append(translation_table[codon])

    for i in range(3):
        translations[i] = "".join(translations[i])
    return translations


def translate(seq):
    """
    Takes a DNA sequence and it's id and returns a dictionary with the 6 different translations
    adds the character '[' as a stop codon. The reason is because [ is the character after Z

    :param seq: dna sequence
    :return to_return: a dictionary of the 6 different translations
    """

    if not seq:
        return dict()

    # reverse_seq = reverse_complement(seq)
    reading_frames = inner_loop(seq)
    # print(reading_frames)
    # reading_frames_reverse = inner_loop(reverse_seq)
    reading_frames = {1: reading_frames[0], 2: reading_frames[1], 3: reading_frames[2]}

    return reading_frames


def translate_no_stop(seq, seq_name):
    """
    Takes a DNA sequence and it's id and returns a dictionary with the 6 different translations
    but it stops once it hits a stop codon

    :param seq: dna sequence
    :param seq_name: the dna sequence id
    :return to_return: a dictionary of the 6 different translations
    """
    cdef int n, i, j

    if not seq:
        return []

    # reverse complement of a seq
    reverse_seq = reverse_complement(seq)

    # I'm checking all 6 different reading frames (3 forward and 3 reverse)
    reading_frames = [""] * 3
    reading_frames_reverse = [""] * 3

    n = len(seq)
    stop_translation = [False, False, False]
    # todo instead of stopping at the stop codon, simply keep it then split
    #   the names here should be unique and indicate which reading frame
    #   I can then simply take that into account and this info as an extra tag
    for i in range(0, n - n % 3, 3):
        # a simple way to check if I hit a stop codon in one of the frames
        # so I stop that and continue the other frames

        for j in range(3):
            if stop_translation[j]:  # if we hit a stop codon, I don't add amino acids to that reading frame anymore
                continue

            # the codon interval in the sequence
            start = i + j
            end = i + j + 3
            if end > n:  # end of sequence
                continue
            codon = seq[start:end].upper()

            if "N" in codon:
                reading_frames[j] += "X"  # we don't know what amino acid
            elif translation_table[codon] == "_":  # stop codon
                stop_translation[j] = True  # in the next iteration, I skip this frame
                continue
            elif codon not in translation_table:
                print("Warning! Codon {} not found in the translation table".format(codon))
            else:
                reading_frames[j] += translation_table[codon]

    stop_translation = [False, False, False]
    for i in range(0, n - n % 3, 3):
        # doing the same for the reverse
        # probably there's a smarter way than copying all the above
        # I'll think about it later

        for j in range(3):
            if stop_translation[j]:
                continue

            start = i + j
            end = i + j + 3
            if end > n:
                continue
            codon = reverse_seq[start:end].upper()
            if "N" in codon:
                reading_frames_reverse[j] += "X"  # we don't know what amino acid
            elif translation_table[codon] == "_":  # stop codon
                stop_translation[j] = True  # in the next iteration, I skip this frame
                continue
            elif codon not in translation_table:
                print("Warning! Codon {} not found in the translation table".format(codon))
            else:
                reading_frames_reverse[j] += translation_table[codon]

    to_return = dict()
    counter = 0
    for r in reading_frames:
        to_return[seq_name + "_reading_frame_" + str(counter)] = r
        counter += 1

    counter = 0
    for r in reading_frames_reverse:
        to_return[seq_name + "_reverse_reading_frame_" + str(counter)] = r
        counter += 1

    return to_return
