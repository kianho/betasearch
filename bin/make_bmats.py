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

# Regex to match the basename of a valid PDB file.
PDB_BASENAME_RE = \
    re.compile(r"(pdb)?([0-9A-Za-z]{4})\.(ent|pdb)(\.gz$)?")


def get_gunzipped_fn(fn):
    """Generated temporary path to a gunzipped copy of a file. The original file
    path is returned if the file is not gzipped.

    Parameters
    ----------
    fn : str
        Path to the gzip file (may not be gzipped).

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
    """Generate the beta-matrices from a given PDB file.

    Parameters
    ----------
    pdb_fn : str
        Path to a PDB file.

    Yields
    ------
    betapy.BetaMatrix
        A beta-matrix object.

    """

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


def write_bmats_fn(bmats_fn, paths_fn=None):
    """Generate beta-matrices from PDB files and write them to a file.

    Parameters
    ----------
    bmats_fn : str
        Path to the file in which the beta-matrices will be written.
    paths_fn : str, optional
        File containing PDB file paths, one per line. If None, the paths will be
        read from stdin (default: None).

    Returns
    -------
    None

    """

    #
    # Read pdb file paths one line at a time from stdin.
    #
    with open(bmats_fn, "wb") as fout:
        if paths_fn:
            fin = open(paths_fn, "rb")
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

    return


if __name__ == "__main__":
    opts = docopt(__doc__)
    write_bmats_fn(opts["--output"], opts["--input"])
