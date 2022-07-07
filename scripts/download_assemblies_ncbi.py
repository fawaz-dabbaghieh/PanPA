import sys, os, time
from Bio import Entrez

'''
This can also be done with a simple bash script, we can loop through
accession numbers and modify this link and download with wget

if in bash our var is the variable for the accession number, this command can be used
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$var&rettype=fasta&retmode=text" -O $var.fasta
'''
Entrez.email = "fawaz@hhu.de"

if len(sys.argv) < 4:
    print("you need to give a list of ncbi accession numbers and output directory")
    sys.exit()

accessions_numbers = []
with open(sys.argv[1], "r") as in_file:
    for l in in_file:
        accessions_numbers.append(l.strip())

for acc in accessions_numbers:
    sequence = []
    handle = Entrez.fetch(db="nucleotide", id=acc,
        rettype="fasta", retmode="text")

    try:
        for l in handles:
            sequence.append(l.strip())
    except:
        pass

    if sequence:
        with open(f"{acc}.fasta", "w") as out_file:
            for s in sequence:
                if s:
                    out_file.write(s + "\n")

    time.sleep(2)  # as too many request per second can cause my ip to be blocked