#!/usr/bin/env python2.7
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    ...

Usage:
    betasearch-local.py -d DIR [-q FILE] [options]

Options:
    -q, --queries FILE              Read one-or-more queries from FILE.
                                    Otherwise read the queries from /dev/stdin.
    -d, --index-dir DIR             Directory containing the betasearch index.
    -o, --output FILE               Write the matching beta-matrices to FILE.
                                    Otherwise write to /dev/stdout by default.
    --validate                      TBD
    --DEBUG                         TBD

"""

import os
import sys
import shelve
import _betasearch as bs
import networkx as nx
import argparse
import pprint
import numpy as np
import shelve

from docopt import docopt
from sys import stderr
from collections import defaultdict
from _betasearch import *

from whoosh.index import open_dir
from whoosh.qparser import QueryParser
from whoosh.query import *

SEP_CHAR = "^"


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

    global SEP_CHAR

    return True, SEP_CHAR + bm_line


def run(query_blob, index_dir):
    r"""Run a single query defined on a single line.

    Parameters
    ----------
        query_blob : str
            The query "blob" of the form: "<query_id> <query_str>".
        index_dir : str
            Path to the betasearch index directory.

    Yields
    ------
        query_id : str
            Unique user-specified query identifier.
        bmat_result : str
            One-line csv representation of a matching beta-matrix.

    """

    vals = query_blob.split(SEP_CHAR)
    query_id, query_str = vals[0], (":" + vals[1])
    lines_db = os.path.join(index_dir, "lines.db")
    shelf = shelve.open(lines_db)

    q = bs.Query(query_str)

    # Whoosh query processing operations.
    whoosh_index = open_dir(os.path.join(index_dir, "whoosh"))
    whoosh_reader = whoosh_index.reader()
    whoosh_searcher = whoosh_index.searcher()

    whoosh_query_parser = QueryParser("trimers", whoosh_index.schema)
    whoosh_query = whoosh_query_parser.parse(q.get_whoosh_query_str())
    whoosh_results = whoosh_searcher.search(whoosh_query, limit=None)

    # Get the matching beta-matrices from disk (using shelve).
    results_gen = \
        (shelf[sheet_id] for sheet_id in q.verify(whoosh_results, index_dir))

    for bmat_result in results_gen:
        yield query_id, bmat_result

    shelf.close()        

    return


if __name__ == "__main__":
    opts = docopt(__doc__)

    if opts["--output"]:
        fout = open(opts["--output"], "wb")
    else:
        fout = sys.stdout

    if opts["--queries"]:
        fin = open(opts["--queries"], "rb")
    else:
        fin = sys.stdin

    # Run multiple queries.
    for query_str in (q.strip() for q in fin):
        for query_id, bmat_result in run(query_str, opts["--index-dir"]):
            fout.write("{} {}".format(query_id, bmat_result) + os.linesep)

    fout.close()
    fin.close()
