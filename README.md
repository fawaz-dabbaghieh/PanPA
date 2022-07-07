# PanPA
PanPA is a tool for building panproteome graphs and aligning sequences back to the graphs.

## Usage
So far the tools is still under development.

For installation, you need to have Cython installed and you can call `python3 setup.py install --user`, which should generate a local binary called panpa that can be used.

PanPA has basically 3 main steps (subcommands):

* Building an index from the input MSA files.
* Building a graph from each MSA.
* Aligning sequences to the graphs generated using the index

```
usage: ProteinAligner [-h] [--log_file LOG_FILE] [--log_level LOG_LEVEL] {build_index,build_gfa,align} ...

Protein Graphs Aligner

Subcommands:
  {build_index,build_gfa,align}
                        Available subcommands
    build_index         Building an index from MSAs
    build_gfa           building GFAs out of MSAs
    align               aligning sequences given to graphs

Global Arguments:
  -h, --help            show this help message and exit
  --log_file LOG_FILE   The name/path of the log file. Default: log.log
  --log_level LOG_LEVEL
                        The logging level [DEBUG, INFO, WARNING, ERROR, CRITICAL]. Default: INFO
```

### Building index
With the `build_index` subcommand, using `-h` will give you the following information:

```
usage: ProteinAligner build_index [-h] [-f IN_FILES [IN_FILES ...]] [-l IN_LIST] [-d IN_DIR] [-o OUT_INDEX]
                                  [--seeding_alg SEEDING_ALG] [-k K-MER] [-w WINDOW] [--seed_limit SEED_LIMIT]

optional arguments:
  -h, --help            show this help message and exit
  -f IN_FILES [IN_FILES ...], --fasta_files IN_FILES [IN_FILES ...]
                        Input MSA(s) in fasta format, one or more file space-separated
  -l IN_LIST, --fasta_list IN_LIST
                        a text file with all input MSAs paths each on one new line
  -d IN_DIR, --in_dir IN_DIR
                        Directory path containing one or more amino acid MSA in FASTA format (gzipped allowed)
  -o OUT_INDEX, --out_index OUT_INDEX
                        The output index file name
  --seeding_alg SEEDING_ALG
                        Seeding algorithm. Choices: k_mers, wk_min. Default: k_mers
  -k K-MER, --kmer_size K-MER
                        K-mer size for indexing the sequencing. Default: 5
  -w WINDOW, --window WINDOW
                        Window size when using w,k-minimizers instead of k-mers for indexing. Default:8
  --seed_limit SEED_LIMIT
                        Indicates how many graphs can a seed belong to. Default: 5, give 0 for no limit
```

The idea is that you can give a directory full of MSAs in fasta format `-d`, a txt file with the list of these input MSAs `-l`, or just give the files as arguments `-f`. User can choose the seeding algorithm like `wk_min` or `k_mers` and the respective sizes of k and/or w.

The `--seed_limit` is the limit for how many graphs can one seed belong to, if 0 then all hits are stored.


### Building GFAs
Similar to the previous subcommand, this takes files, a list of files, or a directory of MSAs that will be turned into graphs in GFA format.

This step can be parallelized, so you can choose how many cores.
```
usage: ProteinAligner build_gfa [-h] [-f IN_FILES [IN_FILES ...]] [-l IN_LIST] [-d IN_DIR] [-c CORES] [-o OUT_DIR]

optional arguments:
  -h, --help            show this help message and exit
  -f IN_FILES [IN_FILES ...], --fasta_files IN_FILES [IN_FILES ...]
                        Input MSA(s) in fasta format, one or more file space-separated
  -l IN_LIST, --fasta_list IN_LIST
                        a text file with all input MSAs paths each on one new line
  -d IN_DIR, --in_dir IN_DIR
                        Directory path containing one or more amino acid MSA in FASTA format (gzipped allowed)
  -c CORES, --cores CORES
                        Numbers of cores to use for aligning
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory where the index files and graphs from the MSAs are stored

```

### Aligning
For aligning, the user need to provide the index that was built with the indexing step, and the graphs built from the MSA, similar to previous the graphs are provided either as individual GFA files, a list or a directory full of the GFAs and sequences to align and cores for parallelization. The user can also use what substitution matrix to use, the gap score and a minimum alignment id cutoff to filter the alignments. The `--seed_limit` parameter basically indicates that after taking the seeds from a read and matching it to which graphs we might want to align to, this limits to how many graphs we align, the seed matches are counted and ordered, meaning if the cutoff is 5, then we align to the top 5 graphs where we mostly had hits with.

```
usage: ProteinAligner align [-h] [-g IN_FILES [IN_FILES ...]] [-l IN_LIST] [-d GRAPHS] [--index INDEX] [-r SEQS]
                            [-c CORES] [--sub_matrix SUB_MATRIX] [-o GAF] [--gap_score GAP_SCORE]
                            [--min_id_score MIN_ID_SCORE] [--seed_limit SEED_LIMIT]

optional arguments:
  -h, --help            show this help message and exit
  -g IN_FILES [IN_FILES ...], --gfa_files IN_FILES [IN_FILES ...]
                        Input GFA graphs, one or more file space-separated
  -l IN_LIST, --gfa_list IN_LIST
                        a text file with all input graphs paths each on one new line
  -d GRAPHS, --in_dir GRAPHS
                        Path to directory with GFA files
  --index INDEX         Path to pickled index file generated in the build step
  -r SEQS, --seqs SEQS  The input sequences to align in fasta format
  -c CORES, --cores CORES
                        Numbers of cores to use for aligning
  --sub_matrix SUB_MATRIX
                        Substitution matrix to use for alignment, default: blosum62
  -o GAF, --out_gaf GAF
                        Output alignments file path
  --gap_score GAP_SCORE
                        The gap score to use for the alignment, default: -3
  --min_id_score MIN_ID_SCORE
                        minimum alignment identity score for the alignment to be outputted, [0,1]
  --seed_limit SEED_LIMIT
                        How many graphs can each seed from the query sequence have hits to, default: 3

```