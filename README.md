# BetaSearch

**WARNING**: We are currently in the process of refactoring BetaSearch, as such
some of the documentation is still incomplete.  You may contact me via email
(see `setup.py`) and I will assist you with your queries.  Alternatively, you
can refer to the original betasearch website
([http://betasearch.kianho.net](http://betasearch.kianho.net)) whilst this
repository is being updated. Thank you for your patience!

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
tests/pdb_files/1ubq.pdb
tests/pdb_files/pdb1brp.ent
tests/pdb_files/pdb1mtp.pdb
tests/pdb_files/pdb3cd4.ent.gz
...
```
then the resulting beta-matrices will be written to `BMATS.dat`
in the following single-line format:
```
sheet-3cd4-A-003^4^5^4^f^f^15^t cell surface glycoprotein cd4^human^homo sapiens^KKVEF,QLVTC,.SVQ.,..GQ.
sheet-3cd4-A-002^2^2^2^f^f^4^t cell surface glycoprotein cd4^human^homo sapiens^LL,VL
sheet-3cd4-A-001^3^3^3^f^f^8^t cell surface glycoprotein cd4^human^homo sapiens^VEL,IIL,AD.
sheet-3cd4-A-000^6^18^8^t^f^51^t cell surface glycoprotein cd4^human^homo sapiens^......KVVLGK......,.QKE-EVQLLVFGLTA..,.VEC-IYTD...ELTLTL,.FHW-KN.......TLSV,QNGLIK............,FLT...............
sheet-3cd4-A-004^2^2^2^f^f^4^t cell surface glycoprotein cd4^human^homo sapiens^GT,ID
sheet-1ubq-A-000^5^8^5^f^f^24^ubiquitin^human^homo sapiens^...TITLE,..TKVFIQ,LVLHLT..,QRLIF...,...QK...
sheet-1brp-A-000^16^33^9^t^t^152^retinol binding protein^human^homo sapiens^.....................CAD-SYSF-VFS,.....................L-LRCSYQ-VAY,......................GNDDHWIVDT.,....................GWYKMKFK.....,...............DVCADMVGTFT.......,...............RVRGKATASM........,..............R-L..VAEFSV........,.............KK-AMAYWTG..........,.......CAD-SYSF-V-FS.............,.......L-LRCSYQ-V-AY.............,........GNDDHWIVD-T..............,......GWYKMKFK...................,.DVCADMVGTFT.....................,.RVRGKATASM......................,.RL..VAEFSV......................,KKAMAYWTG........................
sheet-1mtp-A-001^6^15^6^f^f^62^serine proteinase inhibitor (serpin), chain a^N/A^thermobifida fusca^.FELTQPHQ......,DVRLRAQHIVTDV-Y,GAEGAAATAAM-M-L,RAKAWLANTLI-ARL,...LASRTVLW-V..,........DVR-T..
...
...
```
Generate a BetaSearch index (stored in a directory) from beta-matrices:
```
python ./bin/make_index.py -d ./MY-INDEX/ < BMATS.dat
```
Query a BetaSearch index:
```
python ./bin/run_queries.py -d ./MY-INDEX/ -q QUERIES.txt -o MATCHES.txt
```
where `QUERIES.txt` contains one-or-more beta-matrix queries:
```
query-000 TLE,.I.
query-001 TITLE
query-002 SVT..,SYN..,..SFH
```
The matching beta-matrices in `MATCHES.txt` are formatted as:
```
query-000 sheet-1fyt-D-010^5^8^6^t^f^23^t-cell receptor alpha chain^human^homo sapiens^SVT.....,SYN..LLV,..SFHLTK,..KFEAE.,....LVK.
query-000 sheet-1ubq-A-000^5^8^5^f^f^24^ubiquitin^human^homo sapiens^...TITLE,..TKVFIQ,LVLHLT..,QRLIF...,...QK...
query-000 sheet-1fyt-A-001^4^10^5^t^f^29^hla class ii histocompatibility antigen, dr alpha chain^human^homo sapiens^..EVTVLT..,FKDIFCILVN,F-RKFHYLPF,P-L..ES...
query-001 sheet-1ubq-A-000^5^8^5^f^f^24^ubiquitin^human^homo sapiens^...TITLE,..TKVFIQ,LVLHLT..,QRLIF...,...QK...
query-002 sheet-1fyt-D-010^5^8^6^t^f^23^t-cell receptor alpha chain^human^homo sapiens^SVT.....,SYN..LLV,..SFHLTK,..KFEAE.,....LVK.
...
...
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
