import os
import sys
import pdb


def cigar_to_dict(cigar):
        """
        turns a cigar string into a dict with letter as key and number of base pairs as value
        """
        cigar_dict = dict()  
        number = ""
        for c in cigar:
                if c.isnumeric():
                        number += c
                else:  # we finished the counting for that class
                        if c not in cigar_dict:
                                cigar_dict[c] = int(number)
                        else:
                                cigar_dict[c] += int(number)
                        number = ""
        return cigar_dict


if len(sys.argv) < 3:
        print("you need to give the input sam file and the threshol in %")
        sys.exit()

in_file = sys.argv[1]
thresh = int(sys.argv[2])

with open(in_file, "r") as infile:
        for l in infile:
                if l.startswith("@"):
                        print(l.strip())
                else:
                        line = l.strip().split()
                        if line[1] == "4":  # not aligned filtered out
                                continue
                        cigar_dict = cigar_to_dict(line[5])
                        original_seq_len = 0
                        aligned = 0
                        seq_in_sam = len(line[9])
                        if "H" in cigar_dict:  # there's a hard clip, so the original seq len is seq_in_sam + cigar_dict["H"]
                                original_seq_len = seq_in_sam + cigar_dict["H"]
                                if "S" in cigar_dict:  # there's soft clipping that need to be removed
                                        aligned = original_seq_len - cigar_dict["H"] - cigar_dict["S"]
                                else:
                                        aligned = original_seq_len - cigar_dict["H"]
                        else:
                                original_seq_len = seq_in_sam  # nothing clipped from the sequence in sam
                                if "S" in cigar_dict:
                                        aligned = original_seq_len - cigar_dict["S"]
                                else:  # it is all aligned
                                        aligned = original_seq_len
                        if (aligned * 100)/original_seq_len > thresh:
                                print(l.strip())
                        else:
                                pass
