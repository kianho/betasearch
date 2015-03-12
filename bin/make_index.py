#!/usr/bin/env python
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    Build a betasearch index from beta-matrices. The index is written to disk as
    a directory using the Whoosh and shelve libraries. Beta-matrices are read
    from stdin, one line at a time.

Usage:
    make_index.py -d INDEXDIR [-n NUM] [-l MB] [-v]

Options:
    -d, --index-dir INDEXDIR    The location of the index directory to be built.
    -n, --nprocs NUM            The number of processes to use [default: 1].
    -l, --limit-mb MB           The max. memory to use when building the index [default: 256].
    -v, --verbose               Display verbose output to stderr e.g.
                                diagnostic messages.

"""

import os
import sys
import cPickle
import time
import datetime
import numpy
import shelve
import commands

from docopt import docopt
from collections import defaultdict
from whoosh.fields import Schema, TEXT, STORED
from whoosh.index import create_in
from betasearch import trimers_gen, update_disk_record


if __name__ == "__main__":
    opts = docopt(__doc__)

    if not os.path.exists(opts["--index-dir"]):
        os.mkdir(opts["--index-dir"])

    whoosh_dir = os.path.join(opts["--index-dir"], "whoosh")
    trimers_dir = os.path.join(opts["--index-dir"], "trimers")
    lines_db = os.path.join(opts["--index-dir"], "lines.db")

    if not os.path.isdir(whoosh_dir):
        os.mkdir(whoosh_dir)

    if not os.path.isdir(trimers_dir):
        os.mkdir(trimers_dir)

    # Instantiate blank Whoosh index objects.
    schema = Schema(sheet_id=STORED,
            molecule_name=STORED,
            organism_common_name=STORED,
            organism_scientific_name=STORED,
            trimers=TEXT(phrase=False))
    index_ = create_in(whoosh_dir, schema)
    writer = index_.writer(procs=int(opts["--nprocs"]),
                limitmb=float(opts["--limit-mb"]))

    start_time = time.time()

    shelf = shelve.open(lines_db)

    # TODO: add formal description of beta-matrix format...
    # Read and index each beta-matrix from stdin.
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
            # TODO:
            #   - needs a more scalable solution due to possible filesystem
            #   limits of number of files in a directory, use no-SQL solution
            #   perhaps?
            record_fn = os.path.join(trimers_dir, "%s.record" % sheet_id)
            with open(record_fn, 'w') as f_rec:
                cPickle.dump(record, f_rec)

        except Exception, e:
            sys.stderr.write("ERROR: %s%s" % (e, os.linesep))

    writer.commit()
    shelf.close()

    finish_time = time.time()

    if opts["--verbose"]:
        elapsed_delta = datetime.timedelta(seconds=finish_time - start_time)

        # TODO: use os stat module instead
        approx_size = commands.getoutput("du -hs %s" % opts["--index-dir"]).split()[0]

        sys.stderr.write("elapsed time: %r%s" % (elapsed_delta, os.linesep))
        sys.stderr.write("approx. index size: %s%s" % (approx_size, os.linesep))
