import sys
import os
import matplotlib.pyplot as plt
import upsetplot
import pickle
from math import log


#################### table columns idx constants
A_BWA = 0
NUMB_BWA = 1
LEN_BWA = 2
ID_BWA = 3

A_GA = 4
NUMB_GA = 5
LEN_GA = 6
ID_GA = 7

A_PANPA = 8
NUMB_PANPA = 9
LEN_PANPA = 10
ID_PANPA = 11
####################

if len(sys.argv) < 3:
	print("You need to give the input table pickle and a prefix for the plot files")
	sys.exit()

in_table = sys.argv[1]
plots_prefix = sys.argv[2]

with open(in_table, "rb") as infile:
	big_table = pickle.load(infile)

del big_table['header']
alignment_counts = [0] * 8
combinations = [[], ['BWA'], ["GraphAligner"], ["PanPA"], ["BWA", "GraphAligner"],
["BWA", "PanPA"], ["GraphAligner", "PanPA"], ["BWA", "GraphAligner", "PanPA"]]

id_dist = [[], [], []]
len_dist = [[], [], []]

for seq, info in big_table.items():
	if (info[A_BWA], info[A_GA], info[A_PANPA]) == (0,0,0):  # no alignments
		alignment_counts[0] += 1
	elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (1,0,0):  # only BWA
		alignment_counts[1] += 1
	elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (0,1,0):  # only GA
		alignment_counts[2] += 1
	elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (0,0,1):  # only PanPA
		alignment_counts[3] += 1
	elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (1,1,0):  # BWA - GA
		alignment_counts[4] += 1
	elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (1,0,1):  # BWA - PanPA
		alignment_counts[5] += 1
	elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (0,1,1):  # BA - PanPA
		alignment_counts[6] += 1
	elif (info[A_BWA], info[A_GA], info[A_PANPA]) == (1,1,1):  # all of them
		alignment_counts[7] += 1
		id_dist[0].append(info[ID_BWA])
		id_dist[1].append(info[ID_GA])
		id_dist[2].append(info[ID_PANPA])
		len_dist[0].append(info[LEN_BWA])
		len_dist[1].append(info[LEN_GA])
		len_dist[2].append(info[LEN_PANPA])


tf_frame = upsetplot.from_memberships(combinations, data=alignment_counts)
upsetplot.plot(tf_frame, sort_by="cardinality")
current_figure = plt.gcf()

plt.title("Alignments comparison", loc='left', fontsize="15")
current_figure.savefig(plots_prefix + "_upset_plot.svg", dpi=600)
current_figure.savefig(plots_prefix + "_upset_plot.png", dpi=600)

plt.clf()

# plotting alignment id distribution
# plt.figure(figsize=(8,6))
bins1 = int(1 + (3.22 * log(len(id_dist[0]))))
bins2 = int(1 + (3.22 * log(len(id_dist[1]))))
bins3 = int(1 + (3.22 * log(len(id_dist[2]))))
plt.hist(id_dist[0], bins=bins1, alpha=0.5, label="BWA (DNA)")
plt.hist(id_dist[1], bins=bins2, alpha=0.5, label="GraphAligner (DNA)")
plt.hist(id_dist[2], bins=bins3, alpha=0.5, label="PanPA (AA)")
plt.title("Alignments identity score distributions")
plt.xlabel("Alignment identity score")
plt.ylabel("Frequency")
plt.legend(loc='upper right')
plt.savefig(plots_prefix + "_alignments_id_distribution.svg", dpi=600)
plt.savefig(plots_prefix + "_alignments_id_distribution.png", dpi=600)

plt.clf()

# plotting alignment length percentage distribution
# plt.figure(figsize=(8,6))
bins1 = int(1 + (3.22 * log(len(len_dist[0]))))
bins2 = int(1 + (3.22 * log(len(len_dist[1]))))
bins3 = int(1 + (3.22 * log(len(len_dist[2]))))
plt.hist(len_dist[0], bins=bins1, alpha=0.5, label="BWA (DNA)")
plt.hist(len_dist[1], bins=bins2, alpha=0.5, label="GraphAligner (DNA)")
plt.hist(len_dist[2], bins=bins3, alpha=0.5, label="PanPA (AA)")
plt.title("Percentage alignmened length distributions")
plt.xlabel("Aligned length percentage")
plt.ylabel("Frequency")
plt.legend(loc='upper right')
plt.savefig(plots_prefix + "_alignments_len_percentage_distribution.svg", dpi=600)
plt.savefig(plots_prefix + "_alignments_len_percentage_distribution.png", dpi=600)


plt.clf()
