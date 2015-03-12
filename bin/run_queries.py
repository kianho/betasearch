#!/usr/bin/env python2.7
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    ...

Usage:
    run_queries.py -d DIR [-q FILE] [options]

Options:
    -q, --queries FILE              Read one-or-more queries from FILE.
                                    Otherwise read the queries from /dev/stdin.
    -d, --index-dir DIR             Directory containing the betasearch index.
    -o, --output FILE               Write the matching beta-matrices to FILE.
                                    Otherwise write to /dev/stdout by default.
    -H, --human-readable            Human-readable results.
    --DEBUG                         TBD

"""

import os
import sys
import shelve
import betasearch

from docopt import docopt

from whoosh.index import open_dir
from whoosh.qparser import QueryParser
from whoosh.query import *

SEP_CHAR = "^"


def validate_query(query_str):
    r"""Validate a query string, ensuring that it is well-formed.

    Parameters
    ----------
        query_str : str

    Returns
    -------
        ...

    """

    if len(query_str) == 0:
        return False, "Error: the query cannot be empty"

    bm_lines = [ row.rstrip().replace(" ", '.') for row in query_str.split(",") ]

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

    if "*" in query_str:
        if re.subn(r"\*", "~", query_str)[1] > 2:
            return False, "Error: the query cannot contain more than two wildcard '*' characters."

    bm_line = ",".join(bm_lines)
    g, mat = betasearch.make_graph(bm_line)

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
                if betasearch.get_neighbour_col(row, col + 1, 1, mat) == None:
                    return False, "Error: the query cannot contain a hanging '-' character."

    # 7. check if there are atleast three residues.
    if g.number_of_nodes() < 3:
        return False, "Error: the query needs to contain at least three connected residues."

    # 8. check if the implied graph consists of a single connected component.
    if not betasearch.graph_is_connected(g):
        return False, "Error: all the trimers in the query need to overlap."

    return True, bm_line


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

    vals = query_blob.split()
    query_id, query_str = vals[0], vals[1]
    is_valid, query_str = validate_query(query_str)

    lines_db = os.path.join(index_dir, "lines.db")
    shelf = shelve.open(lines_db)

    q = betasearch.Query(query_str)

    # Whoosh query processing operations.
    whoosh_index = open_dir(os.path.join(index_dir, "whoosh"))
    whoosh_reader = whoosh_index.reader()
    whoosh_searcher = whoosh_index.searcher()

    whoosh_query_parser = QueryParser("trimers", whoosh_index.schema)
    whoosh_query = whoosh_query_parser.parse(q.get_whoosh_query_str())
    whoosh_results = whoosh_searcher.search(whoosh_query, limit=None)

    print list(q.verify(whoosh_results, os.path.join(index_dir, "trimers")))

    # Get the matching beta-matrices from disk (using shelve).
    results_gen = \
        (shelf[metadata[0]] for metadata in
                q.verify(whoosh_results, os.path.join(index_dir, "trimers")))

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
            if opts["--human-readable"]:
                cols = bmat_result.split(SEP_CHAR)
                meta_data = SEP_CHAR.join(cols[:-1])
                bmat_rows = os.linesep.join(cols[-1].split(","))
                output_blob = meta_data + os.linesep + bmat_rows + os.linesep * 2
            else:
                output_blob = query_id + " " + bmat_result + os.linesep

            fout.write(output_blob)

    fout.close()
    fin.close()
