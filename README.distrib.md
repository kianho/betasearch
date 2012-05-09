# README

This .tar.gz file contains the source code for the Python implementation of
BetaSearch, which can be used to build and query your own beta-sheet indices.  

## Help

Feel free to email me
([`hohkhkh1@csse.unimelb.edu.au`](mailto:hohkhkh1@csse.unimelb.edu.au)) should
you have problems running this code or if you have any questions about our work.

## Web server

Alternatively, you can use BetaSearch via our web-server
([here](http://betasearch.servehttp.com/query)).

## Instructions

### Dependencies

* **Operating System**:

    * Linux/Unix (BetaSearch was developed and tested on Ubuntu 11.10 and
      12.04)  

* **Software**:

    * [Python 2.7.x](http://python.org/download/releases/2.7.3/)

                $ sudo apt-get install python2.7

    * [python-networkx](http://networkx.lanl.gov/)

                $ sudo apt-get install python-networkx

    * [python-whoosh](https://bitbucket.org/mchaput/whoosh/wiki/Home)

                $ sudo apt-get install python-whoosh

    * [python-numpy](http://numpy.scipy.org/)

                $ sudo apt-get install python-numpy


### Installation

1. Unpack the `.tar.gz` file:

        $ tar -zxf betasearch-python.tar.gz

   which extract the top directory (`betasearch-py-local`) into the current directory.
   
2. Change into this directory in order to perform the remaining tasks in this
   section:

        $ cd ./betasearch-py-local


### Building a BetaSearch index

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

**NOTE**: Building a BetaSearch index from an entire PDB repository (e.g. ~75K
PDB structures) will take a couple of hours to complete. We have therefore made
available the following precomputed `.bmats` files for download:

* [`pdb2011.bmats.gz`](http://www.csse.unimelb.edu.au/~hohkhkh1/betasearch/files/pdb2011.bmats.gz)

    - generated from the January 3, 2011 snapshot of the PDB.

* [`pdb2012.bmats.gz`](http://www.csse.unimelb.edu.au/~hohkhkh1/betasearch/files/pdb2012.bmats.gz)

    - generated from the February 15, 2012 snapshot of the PDB.

* [`astral95.bmats.gz`](http://www.csse.unimelb.edu.au/~hohkhkh1/betasearch/files/astral95.bmats.gz)

    - generated from the [ASTRAL SCOP 95% sequence ID filtered
      subset](http://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-sel-gs-bib-95-1.75A.tgz).


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
