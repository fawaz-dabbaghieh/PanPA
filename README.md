- [PanPA](#panpa)
  * [Introduction](#introduction)
  * [Installation](#installaion)
  * [Subcommands](#subcommands)
    + [Building Index](#building-index)
    + [Building Graphs](#building-graphs)
    + [Aligning](#aligning)

# PanPA

PanPA (Pan-Proteome Aligner) is a tool written in Cython that builds protein graphs from MSAs, builds an index for the MSAs, and aligns query sequences back to the graphs generated. It is designed to work on amino acid graphs and the alignment can be done using many possible substitution matrix instead of only doing alignment using edit distance. It can also align DNA sequences back to amino acid graphs by translating the DNA sequences into 6 different possible reading frames.

## Introduction
PanPA takes as input any number of MSAs in FASTA format, where each MSA represents one protein, one protein cluster, or one protein family. PanPA has three main steps:

* Indexing MSAs with `build_index` subcommand
* Building Graphs from MSAs in GFA format with `build_gfa` subcommand
* Aligning a query sequence back to the graphs using the `align` subcommand.

In the following sections, installation will be explained and how to use each subcommand

## Installation
PanPA is written in Cython and requires Python 3.6 or higher, and requires Cython.

You can install Cython through pip with and has been tested with this version of Cython `pip install Cython==0.29.21`

Once the requirements are there, then PanPA can be simply installed in the system with `python3 setup.py install --user`

## Subcommands
Both `build_index` and `build_gfa` take any number of MSAs as input, where the MSAs can be given in a directory, a text
with a list of FASTA files paths, or the file paths can be given in the command argument. Same
goes for the `align` subcommand, where it takes the input graphs also in a directory, a list, or individually in the command


### Building Index
The subommands `build_index` takes MSAs as input, and for each sequence in the MSA, seeds are extracted, the user can specify two types of seeds
k-mers or (w,k)-minimizers using the argument `--seeding_alg` which takes either `kmers` or `wk_min`, then
the user needs to specify the k size with `-k, --kmer_size` and w size with `-w, --window`.
The user also needs to give an output file name/location.

The `--seed_limit` argument takes an integer, which specifies a limit to how many MSAs can one seed belong to.
E.g. one k-mer can be present in all MSAs given, the user can specify a limit on that, and the mathces are ordered
based on how many times that seed was present in that MSA and the top `n` will be taken. If the user chooses to keep
all hits, then `0` is given to this argument and all seed hits will be kept in the index.
```
usage: PanPA build_index [-h] [-f IN_FILES [IN_FILES ...]] [-l IN_LIST] [-d IN_DIR] [-o OUT_INDEX]
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

### Building Graphs
Same as `build_index`, `build_gfa` takes the same set of MSAs and an output directory, and for each MSA it builds a corresponding
graph in GFA format, with path lines `p` corresponding to each sequence in the MSA. The 
graphs will be named the same as the original MSAs only with the extension changed to .gfa

This subcommand can be sped up by given it more than one core with `-c, --cores`

NOTE: Once the GFA files are produced, you shouldn't change their names, as they need to be matched correctly
in the index.

```
usage: PanPA build_gfa [-h] [-f IN_FILES [IN_FILES ...]] [-l IN_LIST] [-d IN_DIR] [-c CORES] [-o OUT_DIR]

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
For aligning query sequences to the graphs, you need to give three main inputs to the subcommand `align`:
the index that was built with `--index`, the input graphs which can be a directory, a text file with list, or
given directly in the command, and finally the query sequences in FASTA. If DNA sequences
are given, then the user needs to use the flag `--dna`. The user can also
specify the substitution matrix to use for the alignment, or print a list of possible matrices with
`--sub_matrix_list`. The user can also specify a certain gap score with `--gap_score`, a cutoff on alignment id with
`--min_id_score`, and can set a limit to how many graphs to align to with `--seed_limit`.

This step can be made faster by giving more cores.

The output alignment are in GAF format. To learn more about this format please check here
```
usage: PanPA align [-h] [-g IN_FILES [IN_FILES ...]] [-l IN_LIST] [-d GRAPHS] [--index INDEX] [-r SEQS] [--dna]
                            [-c CORES] [--sub_matrix SUB_MATRIX] [--sub_matrix_list] [-o GAF] [--gap_score GAP_SCORE]
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
  --dna                 Give this flag if the query sequences are DNA and not AA
  -c CORES, --cores CORES
                        Numbers of cores to use for aligning
  --sub_matrix SUB_MATRIX
                        Substitution matrix to use for alignment, default: blosum62
  --sub_matrix_list     When given, a list of possible substitution matrices will be given
  -o GAF, --out_gaf GAF
                        Output alignments file path
  --gap_score GAP_SCORE
                        The gap score to use for the alignment, default: -3
  --min_id_score MIN_ID_SCORE
                        minimum alignment identity score for the alignment to be outputted, [0,1]
  --seed_limit SEED_LIMIT
                        How many graphs can each seed from the query sequence have hits to, default: 3

```