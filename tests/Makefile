SHELL:=/bin/bash
BIN:=../bin

./test-index/.ds: ./test.bmats
	mkdir -p $@ && touch $@
	$(BIN)/make_index.py -d $(@D) < $<

./test.bmats:
	find ./pdb_files/ | grep ".\(ent\|pdb\)\(.gz\)\?$$" | $(BIN)/make_bmats.py -o $@

clean:
	rm -rf ./test.bmats ./test-index
