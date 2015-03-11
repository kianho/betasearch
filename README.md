# BetaSearch

Please refer to the original betasearch website
([http://betasearch.kianho.net](http://betasearch.kianho.net)) whilst this
repository is being updated.

## Installation

1. Install the `setuptools` package using `ez_setup.py`:
```
python ./ez_setup.py
```
2. Install the dependencies using `pip`: 
```
pip install numpy biopython networkx Whoosh docopt
```

3. Install the `betasearch` package:
```
python ./setup.py install
```

## Example usage

Generate beta-matrices from PDB files:
```
python ./bin/make_bmats.py -i PDB_PATHS.dat -o BMATS.dat
```
where `PDB_PATHS.dat` contains a list of PDB file paths, for example:
```
/path/to/1ubq.pdb
/path/to/pdb1mtp.ent.gz
/path/to/1BRP.ent
/path/to/1erb.pdb.gz
...
then the resulting beta-matrices will be written to `BMATS.dat`
in the following single-line format:
```
```

Generate a BetaSearch index from beta-matrices:
```
python ./bin/make_index.py
```

Query a BetaSearch index:
```
python ./bin/run_queries.py
```


### Dependencies

BetaSearch was developed on python 2.6 and 2.7+ using the following libraries:
- Numpy, BioPython, Whoosh, NetworkX, docopt, and ptgraph2.

#### ptgraph2

The `ptgraph2` package was developed as part of
[Pro-origami](http://munk.csse.unimelb.edu.au/pro-origami/), an application for
the automatic generation of protein topology cartoons. If you use `ptgraph2`,
please cite their original paper:

- A. Stivala, M. Wybrow, A. Wirth, J. Whisstock and P. Stuckey, "Automatic
  generation of protein structure cartoons with Pro-origami", Bioinformatics,
  vol. 27, no. 23, pp. 3315-3316, 2011.
  [[article](http://dx.doi.org/10.1093/bioinformatics/btr575)]

## Feedback and Troubleshooting
If you are having trouble running or modifying BetaSearch, please submit an
issue, pull request, or email me. My email can be found in the `setup.py` file.
