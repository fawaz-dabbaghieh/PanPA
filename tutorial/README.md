# Tutorial on PanPA (full working example)
This tutorial is an example of how to use PanPA starting from collected annotations from NCBI of 
10 E. coli assemblies and ends with graphs to align against.

## Requirements
For this tutorial, you need to install `PanPA`, `mmseqs` for clustering, and some MSA software 
like `clustalo` for example. You can get `mmseqs` using conda, brew, docker, or simply downloading the
precompiled version with `wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz
`

## Data
In this tutorial, we will be using 10 E. coli assemblies/annotations randomly selected from RefSeq. The list of ftp
links are listed in `ftp_links.txt`.

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

### Step 1: Download annotations
To download the annotations using the FTP links from RefSeq
```
$ bash download_proteins.sh ftp_links.txt
```
This will download 10 the proteins FASTA file for each assembly.

### Step 2: Separating into groups
We can use 9 of these assemblies to generate the protein clusters, hence graphs and use the last 1 to align back
to the grpahs generated. Therefore, we can mix all the proteins form 9 of these assemblies to generate the clusters
and leave one out for the alignment.

* Let's keep one of these FASTA files for the alignments later, this one was chosen randomly
```
$ gzip -cd GCF_000006625.1_ASM662v1_protein.faa.gz > GCF_000006625.1_ASM662v1_protein.fasta && rm GCF_000006625.1_ASM662v1_protein.faa.gz
```

* We can now merge all sequences from the other 9 into one FASTA file
```
$ for f in *faa.gz;do gzip -cd $f >> all_proteins.fasta && rm $f;done
```

* You can use `fasta_fastq_statistics.sh` to calculate simple statistics on any FASTA or FASTQ file. 
However, it only accepts files where each sequence is contained in one line. Therefore, we can use this
one-liner to remove the new lines in the sequence
```
$ awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < all_proteins.fasta | sed '/^$/d' > tmp && mv tmp all_proteins.fasta

$ bash scripts/fasta_fastq_statistics.sh all_proteins.fasta
61979 reads
27500039 total length
443.699 average read length
```

### Step 3: Generating clusters with mmseqs
Now that we have all proteins from the 9 assemblies, we can cluster them using `mmseqs`.
The parameters chosen here are just an example, but this of course can be changed.
```
$ ./mmseqs easy-linclust all_proteins.fasta all_proteins_cluster tmp --min-seq-id 0.4
$ rm -r tmp/
```
After running `mmseqs`, we get several outputs, a table with cluster names and sequences in 
the cluster `all_proteins_cluster_cluster.tsv`, a FASTA file with the representetive sequences 
`all_proteins_cluster_rep_seq.fasta`, and a FASTA file with all sequences `all_proteins_cluster_all_seqs.fasta`.

We need each cluster to be in a separate FASTA file, you can then use `scripts/extract_clusters.py`
which a simple Python script that takes a simple txt file with sequences names and the FASTA file with all
sequences and an output directory, and it outputs the sequences of each cluster in a separate FASTA file:
```
$ cut -f1 all_proteins_cluster_cluster.tsv | uniq > cluster_names.txt
$ mkdir clusters
$ python3 scripts/extract_clusters.py cluster_names.txt all_proteins_cluster_all_seqs.fasta clusters/
```

### Step 4: Generating MSAs from clusters
This, of course, can be done using many different MSA tools, for this tutorial we used `clustalo`, 
where we first move all clusters that contain one sequence because there's nothing to do, then we run
`clustalo` on each clusters to generate an MSA.

```
$ python3 /scripts/alignment_validation/move_1seq_file_to_msa.py clusters msas

$ for f `ls -1 clusters/`;do ./clustalo --in clusters/$f > msas/$f;done
```
This will take some time to run as there are many clusters.

### Step 5: MSA to GFA
Now that we have many MSAs, we can use `PanPA` to generate a graph for each MSA.

```
$ mkdir graphs
$ PanPA build_gfa -d msas/ -c 4 -o graphs
```

The `build_index` subcommand can take several cores and run in parallel, here we gave it 4 cores, and finished
converting all clusters to graphs in about 2 minutes on a standard laptop.

### Step 6: Indexing
We need to also index the MSAs where we use the index to guide the alignment to which graphs to align to
as we have a 1 to 1 equivalency between an MSA and a GFA, if a seed points to e.g. MSA1 then we align to GFA1.
The user can choose several parameters for indexing and can increase or decrease the seed size depending
on the data used.

```
$ PanPA build_index -d msas/ --seeding_alg wk_min -k 5 -w 3 --seed_limit 0 -o index_k5_w3_no_limit.index
```
This step takes a bit more tan 1 minute

### Step 7: Aligning
Finally, we have generated graphs and an index, we can give both of these to the `align` subcommand in
`PanPA` and some query sequences to do the alignments.

```
$ PanPA align -d graphs/ --index index_k5_w3_no_limit.index -r GCF_000006625.1_ASM662v1_protein.fasta --min_id_score 0.5 --seed_limit 30 -c 4 -o GCF_000006625.1_aligned.gaf
```

This subcommand can also take several cores which makes alignment faster. For these parameters
the aligment was done in about a minute.
