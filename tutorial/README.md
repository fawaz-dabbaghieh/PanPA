# Tutorial on PanPA (full working example)
We start with 16 E. coli annotations randomly selected from NCBI

| **Accession**   | **FTP Link**                                                                      |
|-----------------|-----------------------------------------------------------------------------------|
| GCF_000002515.2 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/515/GCF_000002515.2_ASM251v1 |
| GCF_000002725.2 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/725/GCF_000002725.2_ASM272v2 |
| GCF_000002985.6 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235 |
| GCF_000005825.2 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/825/GCF_000005825.2_ASM582v2 |
| GCF_000005845.2 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2 |
| GCF_000006605.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/605/GCF_000006605.1_ASM660v1 |
| GCF_000006625.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/625/GCF_000006625.1_ASM662v1 |
| GCF_000006645.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/645/GCF_000006645.1_ASM664v1 |
| GCF_000006725.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/725/GCF_000006725.1_ASM672v1 |
| GCF_000006745.1 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/745/GCF_000006745.1_ASM674v1 |


Steps:
* Download annotations for the 10 samples
* Take the first 9 to build the panproteome from and the last one for aligning back against the panproteome
* Mix all proteins together from the 9 samples
* Run mmseq2 to get protein clusters
* Separate each cluster into a fasta file
* Run mafft on each cluster to generate MSAs
* Turn each MSA into a graph using PanPA
* Index the MSAs using PanPA
* Align the last sample's proteins against the graphs generated using PanPA


You can prepare the proteins using the following commands:
```
# download proteins
$ bash download_proteins.sh ftp_links.txt

# Keeping one for alignments later
$ gzip -cd GCF_000006625.1_ASM662v1_protein.faa.gz > GCF_000006625.1_ASM662v1_protein.fasta && rm GCF_000006625.1_ASM662v1_protein.faa.gz

# Merge all the others together into one file
$ for f in *faa.gz;do gzip -cd $f >> all_proteins.fasta && rm $f;done

# You can calculate some stats about the fasta file produced using scripts/fasta_fastq_statistics.sh
# But you need to first remove new lines inside the sequence, can easily be done with awk
$ awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < all_proteins.fasta | sed '/^$/d' > tmp && mv tmp all_proteins.fasta

# calculate stats
$ bash scripts/fasta_fastq_statistics.sh all_proteins.fasta
61979 reads
27500039 total length
443.699 average read length

# running mmseqs on all_proteins.fasta to get the clusters
$ ./mmseqs easy-linclust all_proteins.fasta all_proteins_cluster tmp --min-seq-id 0.4
$ rm -r tmp/

# we can separate the original sequences into clusters with each cluster into fasta file using
# /scripts/extract_clusters.py which takes a list of cluster names
# the list of cluster names can be extracted from the .tsv file output
# the fasta file output with all sequences produced from mmseqs
# and an output directory
$ cut -f1 all_proteins_cluster_cluster.tsv | uniq > cluster_names.txt
$ mkdir clusters
$ python3 scripts/extract_clusters.py cluster_names.txt all_proteins_cluster_all_seqs.fasta clusters/


# Now we need to turn the clusters to MSAs
# clusters with one sequence inside can be just moved
$ python3 /scripts/alignment_validation/move_1seq_file_to_msa.py clusters msas

# For the lefover clusters, you can use your favorite multiple sequence alignment to turn them into
# MSAs, here I am using clustalo
$ for f `ls -1 clusters/`;do ./clustalo --in clusters/$f > msas/$f;done
# This might take some time as there are around 8000 clusters, but they have few sequences in them
# so shouldn't take too long

# Now that we have a collection of MSAs
# we can turn each MSA into a graph using PanPA's build_gfa
$ mkdir graphs
$ PanPA build_gfa -d msas/ -c 4 -o graphs
# This should take around 2 minutes with 4 cores

# We can also build the index at the same time
# for this experiment, I have chosen to use w,k minimizers and k = 5 and w = 3
# with no limit on how many hits each seed can get
# This should take around 1 minute to run
$ PanPA build_index -d msas/ --seeding_alg wk_min -k 5 -w 3 --seed_limit 0 -o index_k5_w3_no_limit.index

# fianlly, we can align the sequences we extraced to the graph, using the index to find where to align
# This also takes a bit over a minute, however, big part of this goes to loading the graphs into memor 
# Also finding seeds, then alignment 
$ PanPA align -d graphs/ --index index_k5_w3_no_limit.index -r GCF_000006625.1_ASM662v1_protein.fasta --min_id_score 0.5 --seed_limit 30 -c 4 -o GCF_000006625.1_aligned.gaf
```

```
# You can get mmseqs with conda, docker, brew or just downloading the compiled version
# In this tutorial, this version of MMseqs2 was used: 19dce0330fecae4c2c151256fa62a3a458830072
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz
```