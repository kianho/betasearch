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
DISTRIB_TAR_GZ_FILE=./betasearch-python.tar.gz

#
# TEST VARIABLES
#
GEN_BMATS_SCRIPT=../betapy/gen_bmats.py
TEST_PDB_ID=1ubq
TEST_INDEX=./test-index
TEST_QUERY=$(TEST_PDB_ID).bmats
#TEST_PDB_ID=1k4r
PDB_PATH=~/PDB

nothing:

test: test_build_index

test_build_index:
	echo $(TEST_PDB_ID) | $(GEN_BMATS_SCRIPT) -b ./$(TEST_PDB_ID).bmats \
		-p ~/PDB -n ./$(TEST_PDB_ID).natpairs -t ./$(TEST_PDB_ID).topos
	./build-index.py -i $(TEST_INDEX) < $(TEST_PDB_ID).bmats 
	./betasearch-local.py -i $(TEST_INDEX) --stdin < ./$(TEST_QUERY)

distrib: .distrib

.distrib: $(PY_SRC)
	tar -cvzf $(DISTRIB_TAR_GZ_FILE) $(PY_SRC) $(QUERIES) $(BETASEARCH_QUERIES) $(SCRIPTS) $(PTGRAPH2_SRC) $(BMATS)

clean:
	rm -f $(DISTRIB_TAR_GZ_FILE)
	rm -f ./*.{bmats,natpairs,topos}
	rm -rf $(TEST_INDEX)
