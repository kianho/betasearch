#
#
#
#
#
#
#
#
#
#
#
#

PY_SRC=./build-index.py ./_betasearch.py ./betasearch-local.py ./betamatrix.py
QUERIES=./experiments/data/queries
BMATS=./experiments/data/astral95.bmats
BETASEARCH_QUERIES=./experiments/data/betasearch/astral95_bmats_queries.txt
PTGRAPH2_SRC=./experiments/scripts/betapy/ptgraph2/{*.py,README}
SCRIPTS=./experiments/scripts/_gen_bmats.py ./experiments/scripts/betapy/*.py


distrib: .distrib


.distrib: $(PY_SRC)
	tar -cvzf ./betasearch-python.tar.gz $(PY_SRC) $(QUERIES) $(BETASEARCH_QUERIES) $(SCRIPTS) $(PTGRAPH2_SRC) $(BMATS)
