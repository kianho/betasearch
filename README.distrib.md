% README
% Kian Ho <hohkhkh1@csse.unimelb.edu.au>
% May 4, 2012

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

`./experiments/scripts/_gen_bmats.py` generates the beta-matrices from pdb file
paths. TODO


## Installation

### Dependencies
