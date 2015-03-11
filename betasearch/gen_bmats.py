#!/usr/bin/env python2.7
"""

$Id: gen_bmats.py,v 5ae03e735bb2 2012/07/01 05:53:21 hohkhkh1 $

Description:
    Generate the beta-matrices from a list of pdb structures read from stdin.
    Each line in stdin must either be a pdb id or a file path to a pdb file.
    The beta-matrices are then written to stdout in the following format:

Usage:
    TODO
     
"""

import os
import sys
import re
import tempfile
import betapy
import pprint
import glob

from sys import stderr

TMP_DIR = "/tmp"

# Regex to match the basename of a valid PDB file.
PDB_BASENAME_RE = \
    re.compile(r"(pdb)?([0-9A-Za-z]{4})\.(ent|pdb)(\.gz)?")


def is_gzipped(fn):
    return fn.endswith(".gz")


def get_gzipped_fn(fn):
    if is_gzipped(fn):
        return fn

    gz_fn = os.path.join(TMPDIR, fn + ".gz")
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


def parse_options():
    """Parse the command-line options.

    Returns:
       optparse.OptionParser object with the command-line parameter
       values. 

    """

    import argparse

    usage_str = \
        "python2.7 ./gen_bmats.py < <file containing pdb ids or pdb fns>"    

    parser = argparse.ArgumentParser(usage=usage_str)

    parser.add_argument("-b", "--bmatsfn", 
            help="file in which the beta-matrices will be written, one per line.")
    parser.add_argument("-n", "--natpairsfn",
            help="file in which the native strand pairs will be written.")
    parser.add_argument("-t", "--toposfn",
            help="file in which the native beta-sheet topologies will be written.")

    parser.add_argument("--stdout", default=False, action="store_true")
    parser.add_argument("-r", "--showresids", default=False, action="store_true")
    parser.add_argument("-B", "--halfbarrels", default=False, action="store_true")
    parser.add_argument("-j", "--json", action="store_true", default=False)
    parser.add_argument("-T", "--tempdir", default="/tmp")

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    options = parser.parse_args()

    return options


if __name__ == "__main__":
    options = parse_options()

    TMP_DIR = options.tempdir

    if options.stdout:
        bmats_f = sys.stdout
        natpairs_f = sys.stdout
        topos_f = sys.stdout
    else:
        bmats_f = open(options.bmatsfn, "ab")
        natpairs_f = open(options.natpairsfn, "ab")
        topos_f = open(options.toposfn, "ab")

    #
    # Read pdb file paths one line at a time from stdin.
    #

    #for line in sys.stdin:
    for line in ["1UBQ.pdb"]:
        pdb_fn = line.strip()
        delete_pdb_fn = False

        if not os.path.isfile(pdb_fn):
            sys.stderr.write("ERROR: %s was not found, skipping...%s"
                    % (pdb_fn, os.linesep))
            continue

        if is_gzipped(pdb_fn):
            pdb_fn = get_gzipped_fn(pdb_fn)
            delete_pdb_fn = True

        # NOTE: this script assumes that pdb filenames assume the format
        # specified by PDB_BASENAME_RE.
        pdb_basename = os.path.basename(pdb_fn)
        pdb_basename_match = PDB_BASENAME_RE.match(pdb_basename)

        assert(pdb_basename_match)

        pdb_id = pdb_basename_match.group(2)

        if pdb_id.startswith("pdb"):
            pdb_id = pdb_id.replace("pdb", "")

        for mat in gen_beta_matrices(pdb_fn):
            # Sometimes ptgraph2 will incorrectly detect pdbids and set them to
            # "0000". Check for this and set it to the id implied by the file
            # name of the pdb file.
            mat.sheet_id = mat.sheet_id.replace("-0000-", "-" + pdb_id + "-")

            if options.showresids:
                bm_buf = \
                    mat.get_resids_repr(basic=True,
                                        halfbarrels=options.halfbarrels)
            elif options.json:
                bm_buf = mat.get_json_string()
            else:
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

            natpairs_f.write("%s^%s\n" % (mat.sheet_id, ",".join(buf))) 
            topos_f.write("%s^%s\n" % (mat.sheet_id, str(mat.topology)))
            bmats_f.write(bm_buf + "\n")

        if delete_pdb_fn:
            os.remove(pdb_fn)

    bmats_f.close()
    natpairs_f.close()
    topos_f.close()
