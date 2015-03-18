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
import gzip
import shutil
import betasearch.betapy as betapy

from sys import stderr
from docopt import docopt

TMP_DIR = "/tmp"
VERBOSE = False

# Regex to match the basename of a valid PDB file.
PDB_BASENAME_RE = \
    re.compile(r"(pdb)?([0-9A-Za-z]{4})\.(ent|pdb)(\.gz$)?")


def get_gunzipped_fn(fn):
    """Transparently return a temporary path to a gunzipped file. The original
    file path is returned if the file is not gzipped.

    Parameters
    ----------
    fn : str
        Path to the gzip file.

    Returns
    -------
    str
        Path to the temporary gunzipped file.
    bool
        True if the original file was a gzip file, False otherwise.

    """

    if not fn.endswith(".gz"):
        return fn, False

    # temporary gunzipped file location.
    gunzipped_fn = \
        os.path.join(TMP_DIR, os.path.basename(os.path.splitext(fn)[0]))

    # pythonically gunzip the file into a temporary file.
    with open(gunzipped_fn, "wb") as tmp:
        shutil.copyfileobj(gzip.open(fn), tmp)

    assert(os.path.exists(gunzipped_fn))

    return gunzipped_fn, True


def gen_beta_matrices(pdb_fn):
    tmp_fn = None
    protein = None

    protein, tmp_fn = betapy.parse_pdb_fn(pdb_fn)

    if tmp_fn and os.path.isfile(tmp_fn):
        os.remove(tmp_fn)

    # TODO: add diagnostic stderr message.
    if not protein:
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
            pdb_fn, delete_pdb_fn = get_gunzipped_fn(pdb_fn) 
            pdb_id = pdb_basename_match.group(2)

            if pdb_id.startswith("pdb"):
                pdb_id = pdb_id.replace("pdb", "")

            for mat in gen_beta_matrices(pdb_fn):
                # Sometimes ptgraph2 will incorrectly detect pdbids and set them to
                # "0000". Check for this and set it to the id implied by the file
                # name of the pdb file.
                mat.sheet_id = mat.sheet_id.replace("-0000-", "-" + pdb_id + "-")
                bm_buf = mat.get_bmat_string(pdb_fn)
                fout.write(bm_buf + os.linesep)

            if delete_pdb_fn:
                os.remove(pdb_fn)

        fin.close()
