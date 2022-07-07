# PanPA
PanPA is a tool for building panproteome graphs and aligning sequences back to the graphs.

## Usage
So far the tools is still under development.

For installation, you need to have Cython installed and you can call `python3 setup.py install --user`, which shoul generate a local binary called panpa that can be used.

PanPA has basically 3 main steps (subcommands):

* Building an index from the input MSA files.
* Building a graph from each MSA.
* Aligning sequences to the graphs generated using the index

```
usage: panpa [-h] [--log_file LOG_FILE] [--log_level LOG_LEVEL] {build_index,build_gfa,align} ...

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

