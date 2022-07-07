import os
import pdb
import sys
import subprocess

input_dir = sys.argv[1]
fasta_files = []
for f in os.listdir(input_dir):
    if f.endswith("fasta") or f.endswith("fa"):
        fasta_files.append(f)

for f in fasta_files:
    out_file = input_dir + os.sep + f.split(".")[0] + "_msa.fasta"
    process = subprocess.run("/home/fawaz/.local/bin/clustalo --in {} > {}".format(f, out_file),
                             stdout=subprocess.PIPE, shell=True)