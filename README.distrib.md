# README

This .tar.gz file contains the scripts and datasets I used to run the Python
implementation of BetaSearch in the 3D substructure search comparisons. This
version of BetaSearch can be used to build and query your own beta-sheet
indices. The source code for the web-interface is contained in a separate
.tar.gz file.


## File Descriptions

* `./build-index.py` reads beta-matrices from stdin (preferably directly from a
.bmat file) and generates the appropriate betasearch indices within a specified
directory. This directory is then used as an index to the betasearch query
script `betasearch-local.py` (_see below_).

* `./betasearch-local.py` is the main betasearch script which reads in queries
(line-by-line from stdin) of the form: `<query id>:XXX,XXX,XXX` where each `X`
is an amino acid. TODO

`./betapy/gen_bmats.py` generates beta-matrices from pdb file paths.


## Instructions

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


$Id$
