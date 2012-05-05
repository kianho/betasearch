# README

This .tar.gz file contains the source code for the Python implementation of
BetaSearch. This version of BetaSearch can be used to build and query your own
beta-sheet indices.

## Help

Feel free to contact me: `Kian Ho <hohkhkh1@csse.unimelb.edu.au>` should you
have problems running this code or have any questions about our work.

## File Descriptions

* `build-index.py` reads beta-matrices from stdin (preferably directly from a
.bmat file) and generates the appropriate betasearch indices within a specified
directory. This directory is then used as an index to the betasearch query
script `betasearch-local.py` (_see below_).

* `betasearch-local.py` is the main betasearch script which reads in queries
  (line-by-line from stdin) of the form: `<query id>:XXX,XXX,XXX` where each
`X` is an amino acid. TODO

* `_betasearch.py` is a Python module that implements the indexing, querying,
  and verification algorithms used by BetaSearch.

* `betapy/betapy.py` is a Python module containing various PDB and protein
  structure-related classes and functions that are used by BetaSearch.

* `betapy/gen_bmats.py` generates beta-matrices from pdb file paths.

* `betapy/ptgraph2` contains scripts of a Python module originally developed by
  `Alex Stivala <astivala@csse.unimelb.edu.au>` for `Pro-Origami`.
  BetaSearch makes extensive use of this module.


## Instructions

### Installation

Unpack the `.tar.gz` file:

        $ tar -zxf betasearch-python.tar.gz

which will place the top directory (`betasearch-py-local`) into the `pwd`, `cd`
into this directory in order to perform the remaining tasks in this section:

        $ cd ./betasearch-py-local


### Building a BetaSearch index

**NOTE**: Building a BetaSearch index from an entire PDB repository (e.g. ~75K
PDB structures) will take several hours to complete.

A `.bmats` file needs to be created before the BetaSearch index can be built.
The `.bmats` file is created from a list of absolute paths to pdb files in the
form of `/path/to/PDB/ub/pdb1ubq.ent.gz` (where `1ubq` is a lower case pdb id),
this list is passed to the `gen_bmats.py` script via `stdin` and the
beta-matrices are written to `example.bmats`, one beta-matrix per line: 

        $ find /path/to/PDB/ -name *.ent.gz | \  
            ./betapy/gen_bmats.py -b ./example.bmats \
                                  -t ./example.topos \
                                  -n ./example.natpairs \
                                  -p /path/to/PDB > example.bmats

Use the `build-index.py` script to read a `.bmats` file from stdin, this will
generate the BetaSearch index using the Whoosh Python module as well as a few
other auxilliary disk-based indices. All of the associated index files are
generated in a single directory (e.g. `./example-index/`):

        $ ./build-index.py -i ./example-index < ./example.bmats


### Querying a BetaSearch index

The `betasearch-local.py` script is used for querying. Queries are read from
`stdin`, one per line, and _must_ conform to the following format:

        <query-id>^ABC,.C.,DEF

where `<query-id>` can be any alphanumeric character. The "`^`" character is used
as a field separator; the remaining text represent a "flat" representation of a
(rectangular) beta-matrix in which rows are separated by "`,`" comma characters,
which is equivalent to the beta-matrix:

        ABC
        .C.
        DEF

The "`.`" character denotes an empty cell in the beta-matrix, which inherently
matches any or no amino acids (similar, but not equivalent to, a `bash` "`*`"
glob wildcard character).

#### An example query

Queries can be stored in a single file (e.g. `queries.bmats`) as follows:

        query-001^ABC,.C.,DEF
        query-002^DFS,.GG,XXX
        query-003^TRW,ASD,RRE
        .
        .
        .

which are then read into `betasearch-local.py` via `stdin` (you must specify
the `--stdin` command-line option):

        $ ./betasearch-local.py -i ./example-index --stdin < ./queries.bmats

the query hits are written to `stdout` in the following format:

        <sheet-id> <molecule name + organism names>
        <sheet-id> <molecule name + organism names>
        <sheet-id> <molecule name + organism names>
        .
        .
        .

A more human-readable results output can be obtained by specifying the
`--humanreadable` command-line option, which displays the beta-matrix in
addition to the metadata above:

        $ ./betasearch-local.py -i ./example-index --stdin \
                                --humanreadable < ./queries.bmats

        stdout:

        <sheet-id> <molecule name + organism names>
        ABC
        .C.
        DEF

        <sheet-id> <molecule name + organism names>
        XXXDFSXX
        ..X.GGYY
        .XXXXXXX

        .
        .
        .
        


`$Id$`
