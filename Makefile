#
#
#
#
#
#
#
#

QUERIES=./experiments/data/queries
BMATS=./experiments/data/astral95.bmats
BETASEARCH_QUERIES=./experiments/data/betasearch/astral95_bmats_queries.txt
DISTRIB_TAR=./betasearch-python.tar
BETAPY_DIR=../betapy
BETAPY_DISTRIB_TAR=$(BETAPY_DIR)/betapy-distrib.tar
DISTRIB_README=./README.distrib.md
BETAPY_DIR=../betapy
README=README.html

#
# TEST VARIABLES
#
GEN_BMATS_SCRIPT=../betapy/gen_bmats.py
#TEST_PDB_ID=1k4r
TEST_PDB_ID=1ubq
TEST_INDEX=./test-index
TEST_QUERY=$(TEST_PDB_ID).bmats
PDB_PATH=~/PDB

#
# ASTRAL95
#
TEST_ASTRAL95_FAKE_PDB_ID=2329
TEST_ASTRAL95_QUERY=$(TEST_ASTRAL95_FAKE_PDB_ID).bmats
ASTRAL95_FAKE_PDB_PATH=/vlsci/VR0127/kian/bio/PDB/fake-ASTRAL95-pdbstyle-1.75A

ASTRAL95_BMATS=./datasets/astral95.bmats
ASTRAL95_BMATS_TOPOS=$(ASTRAL95_BMATS).topos
ASTRAL95_BMATS_NATPAIRS=$(ASTRAL95_BMATS).natpairs
ASTRAL95_BMATS_INDEX=./indices/astral95_bmats_index
ASTRAL95_BMATS_JOBNAME=astral95_bmats
ASTRAL95_BMATS_INDEX_JOBNAME=a95_bmats_index

#
# PDB2012
#
PDB2012_PDB_PATH=/vlsci/data/PDB/COORDINATE_DATA/SNAPSHOT/20120215/
PDB2012_BMATS=./datasets/pdb2012.bmats
PDB2012_BMATS_TOPOS=$(PDB2012_BMATS).topos
PDB2012_BMATS_NATPAIRS=$(PDB2012_BMATS).natpairs
PDB2012_BMATS_JOBNAME=pdb2012_bmats

#
# PDB2011
#
PDB2011_PDB_PATH=/vlsci/data/PDB/COORDINATE_DATA/SNAPSHOT/20110103/
PDB2011_BMATS=./datasets/pdb2011.bmats
PDB2011_BMATS_TOPOS=$(PDB2011_BMATS).topos
PDB2011_BMATS_NATPAIRS=$(PDB2011_BMATS).natpairs
PDB2011_BMATS_JOBNAME=pdb2011_bmats



phony: readme

manuscript:
	( cd ./manuscript ; make )

test: test_build_index

test_build_index:
	echo $(TEST_PDB_ID) | $(GEN_BMATS_SCRIPT) -b ./$(TEST_PDB_ID).bmats \
		-p ~/PDB -n ./$(TEST_PDB_ID).natpairs -t ./$(TEST_PDB_ID).topos
	./build-index.py -i $(TEST_INDEX) < $(TEST_PDB_ID).bmats 
	./betasearch-local.py -i $(TEST_INDEX) -H --stdin < ./$(TEST_QUERY)

astral95_queries_test:
	# echo $(TEST_ASTRAL95_FAKE_PDB_ID) | $(GEN_BMATS_SCRIPT) -b ./$(TEST_ASTRAL95_FAKE_PDB_ID).bmats \
	# 	-p $(ASTRAL95_FAKE_PDB_PATH) -n ./$(TEST_ASTRAL95_FAKE_PDB_ID).natpairs -t ./$(TEST_ASTRAL95_FAKE_PDB_ID).topos
	#./build-index.py -i $(TEST_INDEX) < $(TEST_ASTRAL95_FAKE_PDB_ID).bmats 
	head $(ASTRAL95_BMATS) | ./betasearch-local.py -i $(ASTRAL95_BMATS_INDEX) -H --stdin

# Build the ASTRAL95 .bmats file.
build_astral95_bmats:
	find $(ASTRAL95_FAKE_PDB_PATH) -name "*.ent.gz" | $(GEN_BMATS_SCRIPT) -b $(ASTRAL95_BMATS) \
		-p $(ASTRAL95_FAKE_PDB_PATH) -n $(ASTRAL95_BMATS_NATPAIRS) -t $(ASTRAL95_BMATS_TOPOS) 2> /dev/null

# Build the ASTRAL95 betasearch index.
build_astral95_bmats_index:
	./build-index.py -i $(ASTRAL95_BMATS_INDEX) --procs 2 -l 4096 < $(ASTRAL95_BMATS)

# Build the PDB2012 .bmats file.
build_pdb2012_bmats:
	find $(PDB2012_PDB_PATH) -name "*.ent.gz" | while read line ; \
		do \
			pdbid=`basename $$line | cut -c 4-7` ; \
			if [[ `grep -c $$pdbid $(PDB2012_BMATS)` -ne 0 ]] ; \
			then \
				echo "$$pdbid already seen, skipping" >&2 ; \
				continue ; \
			fi ; \
			echo $$line | timeout3 -t 600 $(GEN_BMATS_SCRIPT) -b $(PDB2012_BMATS) \
			-p $(PDB2012_PDB_PATH) -n $(PDB2012_BMATS_NATPAIRS) -t $(PDB2012_BMATS_TOPOS) ; \
	done
#	find $(PDB2012_PDB_PATH) -name "*.ent.gz" | $(GEN_BMATS_SCRIPT) -b $(PDB2012_BMATS) \
#		-p $(PDB2012_PDB_PATH) -n $(PDB2012_BMATS_NATPAIRS) -t $(PDB2012_BMATS_TOPOS) 2> /dev/null

# Build the PDB2011 .bmats file.
build_pdb2011_bmats:
	find $(PDB2011_PDB_PATH) -name "*.ent.gz" | while read line ; \
		do \
			pdbid=`basename $$line | cut -c 4-7` ; \
			if [[ `grep -c $$pdbid $(PDB2011_BMATS)` -ne 0 ]] ; \
			then \
				echo "$$pdbid already seen, skipping" >&2 ; \
				continue ; \
			fi ; \
			echo $$line | timeout3 -t 600 $(GEN_BMATS_SCRIPT) -b $(PDB2011_BMATS) \
			-p $(PDB2011_PDB_PATH) -n $(PDB2011_BMATS_NATPAIRS) -t $(PDB2011_BMATS_TOPOS) ; \
	done

pbs_astral95_bmats:
	runpbs -j $(ASTRAL95_BMATS_JOBNAME) --walltime 6:0:0 --pvmem 8gb -c "make build_astral95_bmats"
pbs_astral95_bmats_index:
	runpbs -j $(ASTRAL95_BMATS_INDEX_JOBNAME) --walltime 2:0:0 \
	   	--pvmem 10gb --procs 2 --tpn 2 -c "make build_astral95_bmats_index"
pbs_pdb2012_bmats:
	runpbs -j $(PDB2012_BMATS_JOBNAME) --walltime 240:0:0 --pvmem 8gb -c "make build_pdb2012_bmats"
pbs_pdb2011_bmats:
	runpbs -j $(PDB2011_BMATS_JOBNAME) --walltime 240:0:0 --pvmem 8gb -c "make build_pdb2011_bmats"

distrib:
	cp $(DISTRIB_README) ./README.md
	tar --exclude="*.tar" \
	   	--exclude="*.tar.gz" \
	   	--exclude="*.pyc" \
		--exclude="*.pyo" \
	   	--exclude="experiments" \
		--exclude="indices" \
		--exclude="datasets" \
	   	--exclude="README.experiments.md" \
	   	--exclude=".hg*" \
		--exclude="manuscript" \
	   	--exclude="*.DS_Store" \
	   	--exclude="bsmodule.py" \
	   	-h -cvf $(DISTRIB_TAR) ../betasearch-py-local
	gzip -f -c $(DISTRIB_TAR) > $(DISTRIB_TAR).gz

readme:
	Markdown.pl ./README.distrib.md > ./README.html

clean:
	rm -f $(DISTRIB_TAR) $(DISTRIB_TAR).gz
	rm -f ./*.{bmats,natpairs,topos}
	rm -rf $(TEST_INDEX)
	rm -f $(ASTRAL95_BMATS_JOBNAME).{e,o}*
	rm -f $(PDB2012_BMATS_JOBNAME).{e,o}*
	rm -f $(ASTRAL95_BMATS_INDEX_JOBNAME).{e,o}*
	rm -f runpbs.{e,o}*
	rm -f $(README)
