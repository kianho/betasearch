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
import pprint
import numpy as np

from sys import stderr
from collections import defaultdict
from _betasearch import *

from whoosh.index import open_dir
from whoosh.qparser import QueryParser
from whoosh.query import *

DEFAULT_INDEX = os.path.expanduser("~/workspaces/betasearch-py-local/experiments/data/betasearch/astral95_bmats_db/")


def get_pdb_id(sheet_id):
    return sheet_id.split("-")[1].lower()


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


def run_query(line, index_dir=DEFAULT_INDEX):
    req_paths = { "whoosh-dir" : os.path.join(index_dir, "whoosh"), 
                  "trimers-dir" : os.path.join(index_dir, "trimers") }

    q = bs.Query(line)

    index_ = open_dir(req_paths["whoosh-dir"])
    reader_ = index_.reader()
    searcher_ = index_.searcher()

    parser_ = QueryParser("trimers", index_.schema)
    query_ = parser_.parse(q.get_whoosh_query_str())
    results_ = searcher_.search(query_, limit=None)

    return q.verify(results_, req_paths["trimers-dir"])


def do_query(query_str, index_dir=DEFAULT_INDEX):
    """
    """

    return run_query(":" + query_str, index_dir)


def parse_options():
    """Parse the command-line options.

    Returns:
       optparse.OptionParser object with the command-line parameter
       values. 

    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--queries")
    parser.add_argument("-s", "--singlequery")
    parser.add_argument("-i", "--indexdir")
    parser.add_argument("--validate", default=False, action="store_true")
    parser.add_argument("--stdin", default=False, action="store_true",
            help="Read queries from stdin.")
    parser.add_argument("--DEBUG", default=False, action="store_true")
    options = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(-1)

    return options


def run(line, index_dir=DEFAULT_INDEX, VALIDATE_QUERY=False, DEBUG=False, QUIET=False):
    """Run a single query defined on a single line.

    Arguments:
        line -- 

    Returns:
        the number of hits (sheets matching the query).

    """

    vals = line.strip().split("^")
    query_id = vals[0]
    bmtext = vals[-1]

    if VALIDATE_QUERY:
        is_valid, err_txt = validate_query(bmtext)

        if not is_valid:
            sys.stderr.write("ERROR: the query isn't valid\n")
            sys.stderr.write("ERRLOG: %s\n" % err_txt) 
            return 0

    try:
        for sheet_id, mol_name, common_name, sci_name in do_query(bmtext, index_dir):
            print sheet_id, mol_name, common_name, sci_name
    except:
        if not QUIET:
            print "%s\tEXIT_FAILURE" % query_id
        return 0 

    return


if __name__ == "__main__":
    options = parse_options()

    if options.stdin:
        #
        # Read the queries line-by-line from stdin.
        #
        for line in sys.stdin:
            line = line.strip()
            run(line, index_dir=options.indexdir)

    elif options.singlequery:
        # Run a single query from the --singlequery command line option
        # argument.
        run(options.singlequery, options.indexdir)

    elif options.queries:
        #
        # Run multiple queries from a single file.
        #
        with open(options.queries, "rb") as f:
            for line in f:
                run(line, options.indexdir)
