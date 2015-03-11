#!/usr/bin/env python
"""

$Id$

Description:
    This script constructs an index for use with betasearch.

"""

import os
import sys
import cPickle
import time
import datetime
import numpy
import shelve
import commands

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
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--index_dir", required=True,
                        help="directory in which to store the indices.")
    parser.add_argument("-p", "--procs", default=1, type=int,
                        help="no. of processors to use (default=1).")
    parser.add_argument("-l", "--limitmb", default=128, type=int,
                        help="amount of RAM to use (default=128).")
    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    options = parser.parse_args()

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

    if not os.path.isdir(whoosh_dir):
        os.mkdir(whoosh_dir)

    if not os.path.isdir(trimers_dir):
        os.mkdir(trimers_dir)

    schema = Schema(sheet_id=STORED,
            molecule_name=STORED,
            organism_common_name=STORED,
            organism_scientific_name=STORED,
            trimers=TEXT(phrase=False))
    index_ = create_in(whoosh_dir, schema)
    writer = index_.writer(procs=options.procs, limitmb=options.limitmb)

    start_time = time.time()

    shelf = shelve.open(lines_db)

    # Read and index each beta-matrix from stdin.
    # TODO: description of beta-matrix format...
    for line in sys.stdin:
        line = line.strip()
        bmat_fields = line.split("^")

        assert (len(bmat_fields) == 11), \
                "ERROR: this line contains an insuffient number of fields (11) -- %s" % line

        sheet_id = bmat_fields[0]
        entries = bmat_fields[-1]
        molecule_name = bmat_fields[7]
        organism_common_name = bmat_fields[8]
        organism_scientific_name = bmat_fields[9]
        mat = numpy.array(map(lambda row : list(row), entries.split(",")))

        shelf[sheet_id] = line

        record = { "trimers" : defaultdict(set),
                   "trimers-list" : [],
                   "col-index" : {},
                   "row-index" : {} }

        for trimer in trimers_gen(mat):
            update_disk_record(record, trimer)

        try:
            writer.add_document(sheet_id=u"%s" % sheet_id,
                                molecule_name=u"%s" % molecule_name,
                                organism_common_name=u"%s" % organism_common_name,
                                organism_scientific_name=u"%s" % organism_scientific_name,
                                trimers=u"%s" % " ".join(t.id_str for t in
                                                         record["trimers-list"]))

            f = open(os.path.join(trimers_dir, "%s.record" % sheet_id), 'w')
            cPickle.dump(record, f)
            f.close()

        except Exception, e:
            sys.stderr.write("ERROR: %s\n" % e)

    shelf.close()
    writer.commit()

    finish_time = time.time()

    if options.verbose:
        print
        print "elapsed time:", datetime.timedelta(seconds=finish_time - start_time)
        print "approx. index size:", \
                commands.getoutput("du -hs %s" % options.index_dir).split()[0]
