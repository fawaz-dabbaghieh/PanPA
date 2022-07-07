import os
import sys
from collections import defaultdict
import pdb


def getsubs(loc, s):
    substr = s[loc:]
    i = -1
    while (substr):
        yield substr
        substr = s[loc:i]
        i -= 1


def longestRepetitiveSubstring(r, minocc=3):
    occ = defaultdict(int)
    # tally all occurrences of all substrings
    for i in range(len(r)):
        for sub in getsubs(i, r):
            occ[sub] += 1

    # filter out all substrings with fewer than minocc occurrences
    occ_minocc = [k for k, v in occ.items() if v >= minocc]

    if occ_minocc:
        maxkey = max(occ_minocc, key=len)
        return maxkey, occ[maxkey]
    else:
        raise ValueError("no repetitions of any substring of '%s' with %d or more occurrences" % (r, minocc))


def alignment_match(alignment):
    name = alignment.split("\t")[0]
    name = name.split(".")[0]
    try:
        solution = longestRepetitiveSubstring(name, 2)
    except ValueError:
        solution = ''
    if solution == '':
        return False
    else:
        solution = solution[0]
        solution = solution.split("_")
        if len(solution) >= 4:
            return True


def return_name(alignment):
    # return just the name of the sequence without the graph name attached to it
    # so can be matched later
    name = alignment.split("\t")[0]
    name = name.split(".")[0]
    name = name.split("_")
    for idx, n in enumerate(name):
        if "-" in n:
            break
    if "N" not in name[idx + 1]:
        return "_".join(name[:idx + 2])
    else:
        return "_".join(name[:idx + 1])


match_or_mismatch = sys.argv[2]
if match_or_mismatch not in ["match", "mismatch"]:
    print("You need to choose give the keyword match or mismatch")
    sys.exit()

with open(sys.argv[1], "r") as infile:
    for alignment in infile:
        alignment = alignment.strip()
        # printing the one where the seq did not align to its designated graph
        if match_or_mismatch == "match":
            if alignment_match(alignment):
                print(alignment)
        else:
            if not alignment_match(alignment):
                print(alignment)
