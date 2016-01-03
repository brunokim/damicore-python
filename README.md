# DAMICORE

  DAMICORE is an easy-to-use clustering and classification tool.

  DAMICORE is a pipeline of algorithms which is agnostic to the type of data
thanks to NCD, a compressor-based metric which views any piece of data as
simply a string of bits. This is particularly well-suited for heterogenous
datasets, datasets with difficult characteristics extraction and text datasets.

     Data    Metric      Distance matrix    Simplification Phylogenetic tree
    -----.           .---------------------.              .-----------------.
     x_1 |           |     x_1 x_2 ... x_n |              |   x_1-.  x_n    |
     x_2 |  .-----.  | x_1 d11 d12     d1n |    .----.    |        \  |     |
     x_3 |->| NCD |->| x_2 d21 d22     d2n |--->| NJ |--->|         O-O---- |-->
     ... |  '-----'  | ...                 |    '----'    |        /        |
     x_n |           | x_n dn1 dn2     dnn |              | x_2---'         |
    -----'           '---------------------'              '-----------------'
   
            Community detection        Clusters
    -----.                        .-----------------.
         |                        |  ___       ___  |
         |    .-------------.     | |x_1|     |   | |
    tree |--->| Fast Newman |---->| |x_2| ... |C_m| |
         |    '-------------'     | |x_n|     |   | |
         |                        | '---'     '---' |
    -----'                        '-----------------'


## Dependencies

### Installation
* python 2.7
* igraph library (visit http://igraph.org/c for more informations)
* python-igraph (visit http://igraph.org/python for more informations)

### Execution

  DAMICORE relies on compressors to calculate the distance between a pair of
objects (files).

* gzip (available in almost all \*nix systems)
* bzip2 (available in almost all \*nix systems)
* ppmd (available via packages or at http://ctxmodel.net/files/PPMd)

### Recommended

* FigTree (available in http://tree.bio.ed.ac.uk/software/figtree/)

## Usage

The simplest way to use DAMICORE is to provide as argument a directory
containing all files to be compared:

     $ ./damicore.py examples/texts

This outputs a list of pairs of file name and corresponding cluster index.
For now we are lacking a good tool to visualize this clustering, but there are
other tools that might help. We can output intermediate steps into different
files for analysis:

     $ ./damicore.py examples/texts --ncd-output texts.phylip --format phylip \
     --tree-output texts.newick --graph-image texts.png --output texts.clusters

This outputs the NCD matrix using PHYLIP format, the neighbor-joining tree in
Newick format (readable by FigTree), an image with a graph visualization and
the final clusters into another file.

By default, the script uses PPMd as compressor. If you don't have it installed,
try using gzip or bzip2:

     $ ./damicore.py examples/texts --compressor gzip

For more information on available options, see --help.

## Contact

If you believe you have found a bug, or would like to ask for a feature, please
inform me at brunokim@icmc.usp.br.

## Known TODOs

Lots of things to do! Among them:
* Implement UPGMA as tree joining strategy
* Use other community detection algorithm (Girvan-Newman, for example)
* Create/find a cladogram layout for graph visualization
* Include compressors commonly available for Windows (WinRAR, 7zip, etc.)
* Implement classification tool - receives a training dir and a test dir to
  classify
** Also, automate the sampling process and implement k-fold cross-validation
* Implement verbosity level
* Implement other simplifications strategies that produce a sparse graph from
  the distance matrix instead of a tree (kNN, for example)

## License

This software is licensed under the GPLv2.

Texts examples dataset collected by Francisco José Monaco (monaco@icmc.usp.br)
