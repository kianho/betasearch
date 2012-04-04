#!/usr/bin/env python2.7
"""

Revision: $Id$
Author: Kian Ho <hui.kian.ho@gmail.com>
Created:
Description:
    None

"""

import os
import sys
import shelve
import _betasearch as bs
import networkx as nx
import argparse

from collections import defaultdict
from _betasearch import *

from whoosh.index import open_dir
from whoosh.qparser import QueryParser
from whoosh.query import *

QUERY_LOG = "./queries.log"


def get_pdb_id(sheet_id):
    return sheet_id[:4].lower()


def validate_query(bmtext):
    if len(bmtext) == 0:
        return False, "Error: the query cannot be empty"

    bm_lines = [ row.rstrip().replace(" ", '.') for row in bmtext.split(",") ]

    # find the length of the widest line.
    max_width = len(bm_lines[0])

    for line in bm_lines:
        if len(line) > max_width:
            max_width = len(line)

    # right-pad each line with the required number of "." characters to form a
    # rectangular character matrix.
    bm_lines = [ line.ljust(max_width, ".") for line in bm_lines ]

    # Validate the query beta-matrix.

    # 4. check for valid characters.
    if re.match("^[a-zA-Z\.\-\*]+$", "".join(bm_lines)) == None:
        return False, "Error: the query contains invalid characters, please use alphabetic characters to denote amino acids."

    if "*" in bmtext:
        if re.subn(r"\*", "~", bmtext)[1] > 2:
            return False, "Error: the query cannot contain more than two wildcard '*' characters."

    bm_line = ",".join(bm_lines)
    g, mat = bs.make_graph(bm_line)

    # 5. check to see if the matrix is rectangular.
    if len(mat.shape) != 2:
        return False, "Error: the query needs to be a rectangular matrix."

    # 6. check if all stretched peptides have neighbours.
    num_rows, num_cols = mat.shape

    for row in xrange(num_rows):
        for col in xrange(num_cols):
            if col == 0 and mat[row][col] == '-':
                return False, "Error: the query cannot contain a hanging '-' character."

            if mat[row][col] == '.' or mat[row][col] == '-':
                continue

            if col < num_cols - 1:
                if bs.get_neighbour_col(row, col + 1, 1, mat) == None:
                    return False, "Error: the query cannot contain a hanging '-' character."

    # 7. check if there are atleast three residues.
    if g.number_of_nodes() < 3:
        return False, "Error: the query needs to contain at least three connected residues."

    # 8. check if the implied graph consists of a single connected component.
    if not graph_is_connected(g):
        return False, "Error: all the trimers in the query need to overlap."

    return True, ":" + bm_line


def run_query(line):
    index_type = "normal-indices"
    req_paths = { "whoosh-dir" : os.path.join(index_type, "whoosh"), 
                  "seqs-dir" : os.path.join(index_type, "seqs"),
                  "trimers-dir" : os.path.join(index_type, "trimers"),
                  "lines-db" : os.path.join(index_type, "lines.db"),
                  "chain-info-db" : os.path.join(index_type, "chain-info.db")}

    line = line.upper()

    identical_sheets = defaultdict(list)
    num_sheets = 0
    pdb_ids = set()

    expansions = bs.enforced_wc_expansions_iter(line)
    bm_shelf = shelve.open(req_paths["lines-db"], 'r')
    ch_info = shelve.open(req_paths["chain-info-db"], "r")

    qf = open(QUERY_LOG, "a")

    for exp_line in expansions:
        qf.write(exp_line + "\n")
        q = bs.Query(exp_line)

        index_ = open_dir(req_paths["whoosh-dir"])
        reader_ = index_.reader()
        searcher_ = index_.searcher()

        parser_ = QueryParser("trimers", index_.schema)
        query_ = parser_.parse(q.get_whoosh_query_str())
        results_ = searcher_.search(query_, limit=None)

        for sheet_id in q.verify(results_, req_paths["trimers-dir"]):
            num_sheets += 1
            pdb_ids.add(get_pdb_id(sheet_id))
            chain_info = ch_info[sheet_id[:-2]]
            raw_line = bm_shelf[sheet_id]
            is_barrel = False if raw_line.split(":")[-3] == 'f' else True
            sheet_size = int(raw_line.split(":")[-2])
            sheet_line = raw_line.split(":")[-1]
            identical_sheets[(sheet_size, sheet_line, is_barrel)].append(
                (sheet_id, chain_info["org-name"], chain_info["mol-name"]))

    qf.close()

    bm_shelf.close()
    sorted_results = []

    for key in sorted(identical_sheets.iterkeys()):
        sheet_size, sheet_line, is_barrel = key
        val = identical_sheets[key]
        sorted_results.append((sheet_size, is_barrel, sheet_line, val))

    return sorted_results, num_sheets, len(pdb_ids)


def do_query(query_str):
    """
    """

    is_valid, response_text = validate_query(query_str)
    if is_valid:
        query_str = "\n".join(response_text[1:].split(","))
    else:
        query_str = response_text

    if not is_valid:
        return (query_str, False, ([], 0, 0))

    return (query_str, True, run_query(response_text))


def parse_options():
    """Parse the command-line options.

    Returns:
       optparse.OptionParser object with the command-line parameter
       values. 

    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--query", default=None)
    parser.add_argument("-s", "--simpleoutput", default=False,
        action="store_true")
    options = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(-1)

    return options


def main():
    options = parse_options()

    with open(options.query) as f:
        bmtext = ",".join(x.strip() for x in f.xreadlines())

    query_str, _, (sorted_results, num_sheets, num_pdb_ids) = do_query(bmtext)

    if not options.simpleoutput:
        print "# QUERY RESULTS"
        print "# Query:"
        print "#"

        print os.linesep.join("#\t%s" % l for l in query_str.splitlines())

        print "#"
        print "# No. of sheets: %d" % num_sheets
        print "# No. of pdb ids: %d" % num_pdb_ids
        print 

    for sheet_size, is_barrel, sheet_line, matching_descs in sorted_results:
        if not options.simpleoutput:
            print os.linesep.join(sheet_line.split(","))
            print
            for m in matching_descs:
                print "METADATA" + "," + ",".join(m)

            print
            print "---"
            print


    return


if __name__ == "__main__":
    main()
