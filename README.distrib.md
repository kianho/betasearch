# BetaSearch

## Help

Feel free to email me should you have problems running this code or if you have
any questions about our research.

### Web server

Alternatively, you can use BetaSearch via our [web-server](http://porthos.cis.unimelb.edu.au:2081/).


## Instructions

### Dependencies

* **Operating System**:

    * Linux/Unix (BetaSearch was developed on Ubuntu 11.10 and
      12.04)  

* **Software**:

    **NOTE**: I have included the `check_dependencies.py` script to check if the
    required Python modules have been installed into your `$PYTHONPATH`, it can
    be run by calling `python ./check_dependencies.py` on the command line.

    * [Python 2.7.x](http://python.org/download/releases/2.7.3/)

                $ sudo apt-get install python2.7

    * [python-networkx](http://networkx.lanl.gov/)

                $ sudo apt-get install python-networkx

    * [python-whoosh](https://bitbucket.org/mchaput/whoosh/wiki/Home)

                $ sudo apt-get install python-whoosh

    * [python-numpy](http://numpy.scipy.org/)

                $ sudo apt-get install python-numpy
    
    * [python-biopython](http://www.biopython.org)

                $ sudo apt-get install python-biopython


### Installation

1. Unpack the `.tar.gz` file:

        $ tar -zxf betasearch-src.tar.gz

   which extracts the top directory (`betasearch-py-local`) into the current directory.
   
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

**NOTE**: Building a single `.bmats` file from an entire PDB repository (e.g.
~75K PDB structures) will take a few days to complete. We have therefore made
available precomputed (gzipped) `.bmats` files for download from
[here](http://www.csse.unimelb.edu.au/~hohkhkh1/betasearch).

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
