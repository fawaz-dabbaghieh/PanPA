import os
import sys
import pdb


class Alignment:
	def __init__(self, name, gaf_line, orient):
		self.name = name
		self.gaf_line = gaf_line
		self.genes = list()
		self.start = int(name.split("_")[-2])
		self.end = int(name.split("_")[-1])
		if "-" in orient:
			self.orientation = "-"
		else:
			self.orientation = "+"

	def __eq__(self, other):
		if self.name == other.name:
			return True
		else:
			return False

	def __ne__(self, other):
		if self.name == other.name:
			return False
		else:
			return True

	def __hash__(self):
		return hash(self.name)

	def add_gene(self, gene, start):
		gene = gene.split("_")[2]
		if (gene, start) not in self.genes:
			self.genes.append((gene, start))
			if self.orientation == "+":
				self.genes = sorted(self.genes, key= lambda x: x[1])
			else:
				self.genes = sorted(self.genes, key= lambda x: x[1], reverse=True)
		else:
			pass


if len(sys.argv) < 2:
	print("You need to give an input GAF alignment file")
	sys.exit()


alignments = dict()
# the list will have tuples of gene and alignment location
# l[2] is beginning of alignment location
# l[17].split(":")[-1] is the gene or GFA name
with open(sys.argv[1], 'r') as infile:
	for l in infile:
		l = l.strip().split("\t")
		if l[0] not in alignments:
			alignments[l[0]] = Alignment(l[0], l, l[4])
		alignments[l[0]].add_gene(l[17].split(":")[-1], int(l[2]))

order_al = sorted(list(alignments.values()), key= lambda x: x.start)
# pdb.set_trace()

final_gene_order = []
for a in order_al:
	for g in a.genes:
		if g[0] in final_gene_order:
			continue
		if final_gene_order:
			if final_gene_order[-1] != g[0]:
				final_gene_order.append(g[0])
		else:
			final_gene_order.append(g[0])

print(final_gene_order)
