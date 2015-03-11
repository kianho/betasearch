#!/usr/bin/env python
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    Compute beta-matrices from pdb file paths specified from stdin, one-line-at
    a time. Each beta-matrix is compacted into a single line CSV representation
    for convenience of subsequent processing.

Usage:
    make_bmats.py -o FILE [-i FILE]

Options:
    -h, --help              Show this screen.
    -o, --output FILE       Write the beta-matrices to FILE. 
    -i, --input FILE        Read pdb file paths from FILE, otherwise read file
                            paths from stdin.
     
"""

import os
import sys
import re
import tempfile
import pprint
import glob
import betapy

from sys import stderr
from docopt import docopt


TMP_DIR = "/var/tmp"

# Regex to match the basename of a valid PDB file.
PDB_BASENAME_RE = \
    re.compile(r"(pdb)?([0-9A-Za-z]{4})\.(ent|pdb)(\.gz$)?")


def is_gzipped(fn):
    return fn.endswith(".gz")


def get_gzipped_fn(fn):
    if is_gzipped(fn):
        return fn

    gz_fn = os.path.join(TMPDIR, fn + ".gz")

    # TODO: yuck, must correct this.
    os.system("gzip -f -c %s > %s" % (fn, gz_fn))

    return gz_fn


def gen_beta_matrices(pdb_fn):
    tmp_fn = None
    protein = None

    protein, tmp_fn = betapy.parse_pdb_fn(pdb_fn)

    if tmp_fn and os.path.isfile(tmp_fn):
        os.remove(tmp_fn)

    if protein == None:
        return

    for mat in protein.iter_beta_matrices():
        if mat == None:
            continue

        if len(mat) < 2:
            continue
            
        yield mat

    return


if __name__ == "__main__":
    opts = docopt(__doc__)

    #
    # Read pdb file paths one line at a time from stdin.
    #
    with open(opts["--output"], "wb") as fout:
        if opts["--input"]:
            fin = open(opts["--input"], "rb")
        else:
            fin = sys.stdin
        
        for line in fin:
            pdb_fn = os.path.abspath(os.path.expanduser(line.strip()))

            if not os.path.isfile(pdb_fn):
                sys.stderr.write(
                    "ERROR: %s was not found, skipping...%s" %
                        (pdb_fn, os.linesep))
                continue

            # NOTE: this script assumes that pdb filenames assume the format
            # specified by PDB_BASENAME_RE.
            pdb_basename = os.path.basename(pdb_fn)
            pdb_basename_match = PDB_BASENAME_RE.match(pdb_basename)

            if not pdb_basename_match:
                sys.stderr.write(
                    "ERROR: %s is not a valid PDB file name, skipping...%s" %
                        (pdb_fn, os.linesep))
                continue

            # Remember to delete the temporary gunzipped pdb file if used.
            if is_gzipped(pdb_fn):
                pdb_fn = get_gzipped_fn(pdb_fn)
                delete_pdb_fn = True
            else:
                delete_pdb_fn = False

            pdb_id = pdb_basename_match.group(2)

            if pdb_id.startswith("pdb"):
                pdb_id = pdb_id.replace("pdb", "")

            for mat in gen_beta_matrices(pdb_fn):
                # Sometimes ptgraph2 will incorrectly detect pdbids and set them to
                # "0000". Check for this and set it to the id implied by the file
                # name of the pdb file.
                mat.sheet_id = mat.sheet_id.replace("-0000-", "-" + pdb_id + "-")
                bm_buf = mat.get_bmat_string(pdb_fn)

                # Find the native strand pairs and orientations.
                strand_nums = set()
                pairs = []
                for (s1,s2), is_parallel in mat.orients.iteritems():
                    orient = "p" if is_parallel else "a"
                    pairs.append((s1, s2, orient))
                    strand_nums.add(s1)
                    strand_nums.add(s2)

                # Sometimes the strand numbers for beta-sheets may not be
                # contiguous e.g.  1,2,5,6 instead of 1,2,3,4; this is because its
                # original protein change has multiple beta-sheets.  Make the
                # strand numbers in the native pairs and topologies contiguous.
                # NOTE: the strand numbers will remain unchanged in the python
                # objects themselves.
                contig_strand_nums = {}
                for i, num in enumerate(sorted(strand_nums)):
                    contig_strand_nums[ num ] = i + 1

                # Generate the string representation of the native strand pairs.
                buf = []
                for s1, s2, orient in pairs:
                    contig_s1 = contig_strand_nums[ s1 ]
                    contig_s2 = contig_strand_nums[ s2 ]
                    buf.append("%d:%d:%s" % (contig_s1, contig_s2, orient))

                fout.write(bm_buf + os.linesep)

            if delete_pdb_fn:
                os.remove(pdb_fn)

        fin.close()
