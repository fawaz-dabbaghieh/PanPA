import pickle
import sys

in_table = sys.argv[1]
out_file = sys.argv[2]
with open(in_table, "rb") as in_file:
    big_table = pickle.load(in_file)

header = big_table['header']
del big_table['header']
header = ['sequence_name'] + header

with open(out_file, "w") as outfile:
    outfile.write("\t".join(header) + "\n")
    for s, l in big_table.items():
        string_l = [s] + [str(x) for x in l]
        outfile.write("\t".join(string_l) + "\n")
