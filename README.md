# BetaSearch

## Overview

BetaSearch is a fast method for searching beta-sheets for clusters of
amino acid residues that are connected by peptide- or hydrogen-bonds.

TODO

## Installation

### Dependencies

TODO

## Usage

### Generate beta-matrices from PDB files

The [./bin/make_bmats.py](bin/make_bmats.py) script is used to generate
beta-matrices from PDB files using the following single-line representation,
with fields separated by the `^` character:
```
ls -1 /path/to/pdb/*.{pdb,ent,pdb.gz,ent.gz} | \
                python ./bin/make_bmats.py -o ./bmats.txt
```
output (`bmats.txt`):
```
...
sheet-3cd4-A-003^4^5^4^f^f^15^t cell surface glycoprotein cd4^human^homo sapiens^KKVEF,QLVTC,.SVQ.,..GQ.
sheet-3cd4-A-002^2^2^2^f^f^4^t cell surface glycoprotein cd4^human^homo sapiens^LL,VL
sheet-3cd4-A-001^3^3^3^f^f^8^t cell surface glycoprotein cd4^human^homo sapiens^VEL,IIL,AD.
sheet-3cd4-A-000^6^18^8^t^f^51^t cell surface glycoprotein cd4^human^homo sapiens^......KVVLGK......,.QKE-EVQLLVFGLTA..,.VEC-IYTD...ELTLTL,.FHW-KN.......TLSV,QNGLIK............,FLT...............
sheet-3cd4-A-004^2^2^2^f^f^4^t cell surface glycoprotein cd4^human^homo sapiens^GT,ID
sheet-1ubq-A-000^5^8^5^f^f^24^ubiquitin^human^homo sapiens^...TITLE,..TKVFIQ,LVLHLT..,QRLIF...,...QK...
sheet-1fyt-A-002^4^7^4^f^f^22^hla class ii histocompatibility antigen, dr alpha chain^human^homo sapiens^.LLKHWE,EVRCDYV,NVTWLR.,...VPK.
sheet-1fyt-A-001^4^10^5^t^f^29^hla class ii histocompatibility antigen, dr alpha chain^human^homo sapiens^..EVTVLT..,FKDIFCILVN,F-RKFHYLPF,P-L..ES...
sheet-1brp-A-000^16^33^9^t^t^152^retinol binding protein^human^homo sapiens^.....................CAD-SYSF-VFS,.....................L-LRCSYQ-VAY,......................GNDDHWIVDT.,....................GWYKMKFK.....,...............DVCADMVGTFT.......,...............RVRGKATASM........,..............R-L..VAEFSV........,.............KK-AMAYWTG..........,.......CAD-SYSF-V-FS.............,.......L-LRCSYQ-V-AY.............,........GNDDHWIVD-T..............,......GWYKMKFK...................,.DVCADMVGTFT.....................,.RVRGKATASM......................,.RL..VAEFSV......................,KKAMAYWTG........................
sheet-1mtp-A-001^6^15^6^f^f^62^serine proteinase inhibitor (serpin), chain a^N/A^thermobifida fusca^.FELTQPHQ......,DVRLRAQHIVTDV-Y,GAEGAAATAAM-M-L,RAKAWLANTLI-ARL,...LASRTVLW-V..,........DVR-T..
...
```
where the fields are (in left-to-right order):
- sheet identifier in the format: `sheet-<pdb id>-<chain id>-<ordinal sheet number>`
- number of rows in the beta-matrix
- number of columns in the beta-matrix
- number of beta-strands in the beta-sheet

followed by meta-data obtained from PDB header contents:
- beta-sheet bifurcation: `t` if the beta-sheet is bifurcated, `f` otherwise
- number of residues in the beta-sheet
- molecule name
- organism name (common)
- organism name (scientific)

finally, comma-separated strings denoting the rows of the beta-matrix. For
example, the beta-matrix representation of the beta-matrix in ubiquitin (PDB
ID: [1ubq][1ubq]) is:
```
...TITLE            
..TKVFIQ            
LVLHLT..            
QRLIF...            
...QK...            
```
which has a corresponding single-line representation:
```
...TITLE,..TKVFIQ,LVLHLT..,QRLIF...,...QK...
```

### Generate an index from beta-matrices

The [bin/make_index.py](bin/make_index.py) script is used to create a BetaSearch
index inside as a single directory from a file containing beta-matrices (see
`bmats.txt` output above):
```
python bin/make_index.py -d ./betasearch-index < ./bmats.txt
```

### Query an index for beta-residue motifs

The [bin/run_queries.py](bin/run_queries.py) script is used to run
queries on a BetaSearch index:
```
echo "query-XXX V.,LR" | python bin/run_queries -d ./betasearch-index # single query from stdin
python bin/run_queries -d ./betasearch-index -q ./queries.txt # multiple queries from a file
```
where each query is formatted as a single line:
```
<query id> <one-line beta-matrix representation>
```
for example:
```
query-003 I..,VFI,.T.
```
queries for the motif:
```
I..
VFI
.T.
```
Multiple queries can be performed by appending more single-line queries.


### Citation

Please consider citing our paper if you find BetaSearch useful in your research:

- H. K. Ho, G. Gange, M. J. Kuiper, and K. Ramamohanarao, “BetaSearch: a new method for querying
beta-residue motifs,” _BMC Research Notes_, vol. 5, 2012. [[article][betasearch-doi]]

## Contributors

see [CONTRIBUTORS.md](CONTRIBUTORS.md)

[1ubq]: http://pdb.org/pdb/explore/explore.do?structureId=1ubq
[betasearch-doi]: http://dx.doi.org/10.1186/1756-0500-5-391
