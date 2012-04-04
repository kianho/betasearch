#!/usr/bin/env python
"""

$Id$

Description:
    This script constructs an index for use with betasearch.

"""

import os
import sys
import optparse
import cPickle
import time
import datetime
import numpy
import shelve

from collections import defaultdict
from whoosh.fields import Schema, TEXT, STORED
from whoosh.index import create_in
from _betasearch import trimers_gen, update_disk_record


def parse_options():
    """Parse the command-line options.

    Returns:
       optparse.OptionParser object with the command-line parameter
       values. 

    """

    usage_str = \
        "Usage: %prog"    

    parser = optparse.OptionParser()
    parser.add_option("-i", "--index_dir", default=None,
                        help="directory in which to store the indices.")
    parser.add_option("-p", "--procs", default=1, type="int",
                        help="no. of processors to use (default=1).")
    parser.add_option("-l", "--limitmb", default=128, type="int",
                        help="amount of RAM to use (default=128).")
    parser.add_option("-v", "--verbose", default=False, action="store_true")
    options, _ = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(-1)
    
    return options



if __name__ == "__main__":
    options = parse_options()

    if not os.path.exists(options.index_dir):
        os.mkdir(options.index_dir)

    whoosh_dir = os.path.join(options.index_dir, "whoosh")
    trimers_dir = os.path.join(options.index_dir, "trimers")
    lines_db = os.path.join(options.index_dir, "lines.db")

    os.mkdir(whoosh_dir)
    os.mkdir(trimers_dir)

    schema = Schema(sheet_id=STORED, trimers=TEXT(phrase=False))
    index_ = create_in(whoosh_dir, schema)
    writer = index_.writer(procs=options.procs, limitmb=options.limitmb)

    start_time = time.time()

    shelf = shelve.open(lines_db)

    # Read and index each beta-matrix from stdin.
    # TODO: description of beta-matrix format...
    for line in sys.stdin:
        line = line.strip()
        doc_id, entries = line.split(":")
        mat = numpy.array(map(lambda row : list(row), entries.split(",")))

        shelf[doc_id] = line

        record = { "trimers" : defaultdict(set),
                   "trimers-list" : [],
                   "col-index" : {},
                   "row-index" : {} }

        for trimer in trimers_gen(mat):
            update_disk_record(record, trimer)

        writer.add_document(sheet_id=u"%s" % doc_id,
                            trimers=u"%s" % " ".join(t.id_str for t in
                                                     record["trimers-list"]))

        f = open(os.path.join(trimers_dir, "%s.record" % doc_id), 'w')
        cPickle.dump(record, f)
        f.close()

        if options.verbose:
            print doc_id, "indexed sucessfully"

    shelf.close()
    writer.commit()

    finish_time = time.time()

    print "elapsed time:", datetime.timedelta(seconds=finish_time - start_time)