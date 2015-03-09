#!/usr/bin/env python2.7
"""

$Id: betapy.py,v 9bc7964cd8e5 2012/07/12 15:04:46 hohkhkh1 $

Description:


"""

import os
import sys
import tempfile
import cPickle
import copy
import re
import json
import gzip
import numpy


DSSP_HB_THRESH = -0.5 # -0.5 kcal/mol

try:
    import networkx
except ImportError:
    sys.stderr.write("ERROR: networkx module is not installed\n"
                     "       it can be downloaded from http://networkx.lanl.gov\n")
    sys.exit(1)

try:
    import Bio
except ImportError:
    sys.stderr.write("ERROR: Bio module is not installed\n"
                     "       it can be downloaded from http://www.biopython.org\n")
    sys.exit(1)


from Bio.PDB import calc_dihedral, DSSP
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Vector import Vector
from Bio.SCOP.Raf import protein_letters_3to1

from sys import stderr
from subprocess import Popen, PIPE
from itertools import product, combinations, izip, islice, imap, permutations
from collections import defaultdict, OrderedDict, namedtuple, deque
from ptgraph2 import ptutils, make_graphs


## Monkey patches

@property
def get_bridge_type(self, res_b):
    """

    Returns:
        a tuple of (<"HB"/"nHB">, <"HB"/"nHB">) denoting the
        orientation and whether or not res_a and res_b are hydrogen bonded or
        not.

    """

    def get_ap_hb_pos_types(res_a, res_b):
        """
        """

        if ((res_a.dssp_hb_don1[0] == res_a.dssp_hb_acc1[0]== res_b.dssp_num and
             res_a.dssp_hb_don1[1] <= DSSP_HB_THRESH and 
             res_a.dssp_hb_acc1[1] <= DSSP_HB_THRESH) or 
            (res_a.dssp_hb_don1[0] == res_a.dssp_hb_acc2[0] == res_b.dssp_num and
             res_a.dssp_hb_don1[1] <= DSSP_HB_THRESH and 
             res_a.dssp_hb_acc2[1] <= DSSP_HB_THRESH) or 
            (res_a.dssp_hb_don2[0] == res_a.dssp_hb_acc1[0] == res_b.dssp_num and
             res_a.dssp_hb_don2[1] <= DSSP_HB_THRESH and 
             res_a.dssp_hb_acc1[1] <= DSSP_HB_THRESH) or 
            (res_a.dssp_hb_don2[0] == res_a.dssp_hb_acc2[0] == res_b.dssp_num and
             res_a.dssp_hb_don2[1] <= DSSP_HB_THRESH and 
             res_a.dssp_hb_acc2[1] <= DSSP_HB_THRESH)):
            return ("HB", "HB")

        return ("nHB", "nHB")

    def get_p_hb_pos_types(res_a, res_b):
        """
        """

        if ((res_a.dssp_hb_don1[0] == res_b.dssp_num - 1 and
             res_a.dssp_hb_acc1[0] == res_b.dssp_num + 1 and
             res_a.dssp_hb_don1[1] <= DSSP_HB_THRESH and 
             res_a.dssp_hb_acc1[1] <= DSSP_HB_THRESH) or 
            (res_a.dssp_hb_don1[0] == res_b.dssp_num - 1 and
             res_a.dssp_hb_acc2[0] == res_b.dssp_num + 1 and
             res_a.dssp_hb_don1[1] <= DSSP_HB_THRESH and 
             res_a.dssp_hb_acc2[1] <= DSSP_HB_THRESH) or 
            (res_a.dssp_hb_don2[0] == res_b.dssp_num - 1 and
             res_a.dssp_hb_acc1[0] == res_b.dssp_num + 1 and
             res_a.dssp_hb_don2[1] <= DSSP_HB_THRESH and
             res_a.dssp_hb_acc1[1] <= DSSP_HB_THRESH) or
            (res_a.dssp_hb_don2[0] == res_b.dssp_num - 1 and
             res_a.dssp_hb_acc2[0] == res_b.dssp_num + 1 and
             res_a.dssp_hb_don2[1] <= DSSP_HB_THRESH and
             res_a.dssp_hb_acc2[1] <= DSSP_HB_THRESH)):
            return ("HB", "nHB")

        return ("nHB", "HB")

    return pair_type

@property
def get_pdb_id(self):
    return self.get_full_id()[0].lower()

@property
def get_chain_id(self):
    return self.get_full_id()[2]

@property
def get_pdb_resid(self):
    return ("%s%s" % (self.id[1], self.id[2])).strip()

@property
def get_global_id(self):
    return "%s:%s:%s:%s" % (self.pdb_id, self.chain_id, self.pdb_resid,
                            self.single_aa)

@property
def get_single_aa(self):
    try:
        return three_to_one(self.get_resname())
    except KeyError:
        return protein_letters_3to1(self.get_resname(), "X")

LAST_X1_POINT = defaultdict(
        lambda : "CG",
            { "S" : "OG",
              "T" : "OG1",
              "C" : "SG",
              "I" : "CG1",
              "V" : "CG1" })

@property
def get_X1_angle(self):
    try:
        N = Vector(self["N"].coord)
        CA = Vector(self["CA"].coord)
        CB = Vector(self["CB"].coord)
        last = Vector(self[LAST_X1_POINT[self.single_aa]].coord)

        return (180.0 / numpy.pi) * calc_dihedral(N, CA, CB, last)
    except:
        return None


Bio.PDB.Residue.Residue.pdb_id = get_pdb_id
Bio.PDB.Residue.Residue.chain_id = get_chain_id
Bio.PDB.Residue.Residue.pdb_resid = get_pdb_resid
Bio.PDB.Residue.Residue.global_id = get_global_id
Bio.PDB.Residue.Residue.single_aa = get_single_aa
Bio.PDB.Residue.Residue.X1_angle = get_X1_angle


#
# Code taken from BioPython to simplify residue name lookups.
#

standard_aa_names=["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", 
                   "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",
                   "TRP", "TYR"]

aa1="ACDEFGHIKLMNPQRSTVWY"
aa3=standard_aa_names

d1_to_index={}
dindex_to_1={}
d3_to_index={}
dindex_to_3={}

# Create some lookup tables
for i in range(0, 20):
    n1=aa1[i]
    n3=aa3[i]
    d1_to_index[n1]=i
    dindex_to_1[i]=n1
    d3_to_index[n3]=i
    dindex_to_3[i]=n3


def three_to_one(s):
    """Three letter code to one letter code.
    
    >>> three_to_one('ALA')
    'A'
    >>> three_to_one('TYR')
    'Y'

    For non-standard amino acids, you get a KeyError:

    >>> three_to_one('MSE')
    Traceback (most recent call last):
       ...
    KeyError: 'MSE'
    """

    i=d3_to_index[s]

    return dindex_to_1[i]



#
# Graeme's griddify code
#

# Builds a mapping node -> rank, subject to the matrix.

def build_rowtable(mat):
    """
    """
    rowtable = {}
    for i in range(0,len(mat)):
        for sheet in mat[i]:
            for node in sheet:
                rowtable[node] = i
    return rowtable
    

## Orders the ranks in the sheet
def order_ranks(mat,bridges):
    """
    """

    # Initial row table.
    rowtable = build_rowtable(mat)

    # Generate the edge lists.
    rankedges = [ [] for sheets in mat ]
    for i in range(0,len(mat)):
        dests = set([ rowtable[y] for ((x,y),p) in bridges if rowtable[x] == i])
        for dest in dests:
            if dest not in rankedges[i]:
                rankedges[i].append(dest)
            if i not in rankedges[dest]: # added Nov-3 -- 4:50pm by ggange
                rankedges[dest].append(i)

    # Generate the initial extrema.
    extrema = None
    for i in range(0,len(rankedges)):
        if len(rankedges[i]) == 1:
            extrema = i
            break

    assert extrema != None

    # Construct the order
    ordered = []
    while extrema != None:
       ordered.append(extrema)
       next = rankedges[extrema][0]
       rankedges[next].remove(extrema)

       if len(rankedges[next]) > 0:
           assert len(rankedges[next]) == 1
           extrema = next
       else:
           if next not in ordered: # added Nov-3 -- 4:53pm by kian
               ordered.append(next)
           extrema = None

    assert len(ordered) == len(mat)

    return [ mat[i] for i in ordered ]


## Generates the list of columns that must be aligned.
def gen_cols(nodes,mat,bridges):
    cols = []
    coltable = {}

    nodeset = list(nodes)

    # Order the ranks, and build the new rowtable
    mat = order_ranks(mat,bridges)
    rowtable = build_rowtable(mat)

    ## Generate the set of columns 
    while len(nodeset) > 0:
       ## Fresh column
       node = nodeset.pop()
       
       colid = len(cols)
       cols.append([node])
       coltable[node] = colid
      
       nodeq = [node]
      
       # Follow chains of bonds to generate the rest of the column.
       while len(nodeq) > 0:
          next = nodeq.pop()

          bonded = [ x for ((x,y),p) in bridges if y == next ] + [ y for ((x,y),p) in bridges if x == next ]
          for bond in bonded:
              if bond in nodeset:
                  cols[colid].append(bond)
                  nodeq.append(bond)
                  nodeset.remove(bond)
                  coltable[bond] = colid

    # Sort the columns.  
    for col in cols:
        col.sort(key=lambda node: rowtable[node]) 

    ## Do a topological sort on the columns.
    # Generate the incoming edges.
    inedges = [ 0 for col in cols ]
    outedges = [ [] for col in cols ]
    for sheets in mat:
        sheet = reduce(lambda x, y: x + y, sheets)
        for i in range(0,len(sheet)-1):
            inedges[coltable[sheet[i+1]]] += 1
            outedges[coltable[sheet[i]]].append(coltable[sheet[i+1]])

    # Initialize the set of cols with no incoming edges.
    empty = []
    sorted = []
    for i in range(0,len(inedges)):
        if inedges[i] == 0:
            empty.append(i)

    while len(empty) > 0:
        next = empty.pop()
        sorted.append(cols[next])

        for dest in outedges[next]:
            inedges[dest] -= 1
            if inedges[dest] == 0:
                empty.append(dest)

    assert( len(sorted) == len(cols) )
    return sorted


# Builds a table mapping column number -> grid position.
def build_postable(nodes,mat,bridges):
    """
    """

    cols = gen_cols(nodes,mat,bridges)

    coltable = {}
    for i in range(0,len(cols)):
        for node in cols[i]:
            coltable[node] = i

    # The column of the successor.
    successor = {}
    predecessor = {}
    for sheets in mat:
        sheet = reduce(lambda x, y: x + y, sheets)
        for i in range(0,len(sheet)-1):
            successor[sheet[i]] = coltable[sheet[i+1]]
      
        for onesheet in sheets:
            for i in range(1,len(onesheet)):
                predecessor[onesheet[i]] = coltable[onesheet[i-1]]

    lbtable = {}
    for colid in range(0,len(cols)):
        if colid not in lbtable:
            lbtable[colid] = 0
      
        nextpos = lbtable[colid] + 1

        for node in cols[colid]:
           if node not in successor:
               continue

           if successor[node] not in lbtable or lbtable[successor[node]] < nextpos:
               lbtable[successor[node]] = nextpos

    postable = {}
    for colid in range(len(cols)-1,-1,-1):
        if colid not in postable:
            postable[colid] = lbtable[colid]
      
        nextpos = postable[colid] - 1 
      
        for node in cols[colid]:
            if node not in predecessor:
                continue
            if predecessor[node] not in postable or postable[predecessor[node]] > nextpos:
                postable[predecessor[node]] = nextpos

    return postable


def build_grid(nodes,mat,bridges):
    """
    """

    postable = build_postable(nodes,mat,bridges)
    mat = order_ranks(mat,bridges)
    width = max(postable.values())+1
    cols = gen_cols(nodes,mat,bridges)

    coltable = {}
    for i in range(0,len(cols)):
        for node in cols[i]:
            coltable[node] = i

    grid = []
    for rank in mat:
        row = []
        sheet = reduce(lambda x, y: x + y, rank)
        off = 0

        for i in range(0,width):
            if off < len(sheet) and postable[coltable[sheet[off]]] == i:
                row.append(sheet[off])
                off += 1
            else:
                row.append(None)
      
        grid.append(row)
    return grid


def parse_pdb_fn(pdb_fn, tmp_dir="/tmp", verbose=False):
    tmp_pdb_fn = None

    if os.path.splitext(pdb_fn)[-1] == ".gz":
        pdb_f, tmp_pdb_fn = tempfile.mkstemp(dir=tmp_dir)
        os.system("gunzip -c %s > %s" % (pdb_fn, tmp_pdb_fn))
        os.close(pdb_f)

    if tmp_pdb_fn:
        protein = Protein(tmp_pdb_fn)
    else:
        protein = Protein(pdb_fn)

    return protein, tmp_pdb_fn


#
# CUSTOM EXCEPTIONS
#

class BadSheetMatrixException(Exception):
    """Raised when a sheet matrix is unable to be generated for a sheet.

    """

    pass


class SelfBridgeException(Exception):
    """Raised to signal the occurrence of bridge between two residues within the
    same strand.

    """

    pass


class BridgesNotSatisfiedException(Exception):
    """Raised to signal that not all bridge edges have been satisfied in the raw
    matrix generated by griddify.build_grid. This exception is raised by
    Sheet.build_sheet_matrix(...).

    """

    pass


#
# CLASS DEFINITIONS
#

class DSSPChain(object):
    """A single protein chain parsed from DSSP output.

    """

    def __init__(self, raw_chain):
        """Constructor.

        Args:
            raw_chain: A list of BIO.PDB.Residue objects.

        """

        self.residues = raw_chain
        self.pdb_resid_dict = {}

        for ord_index, res in enumerate(self.residues):
            res.ord_index = ord_index
            self.pdb_resid_dict[res.pdb_resid] = res

        return

    def __iter__(self):
        return iter(self.residues)

    def __getitem__(self, pdb_resid):
        return self.pdb_resid_dict[pdb_resid]

    def get_pdb_id(self):
        return os.path.basename(self.residues[0].full_id[0]).split(".")[0]

    def get_seq(self):
        """Get the one-letter amino acid sequence.

        """
        return "".join(imap(three_to_one,
                            (r.get_resname() for r in self.residues)))

    def get_ss_seq(self):
        """Get the (8-state) DSSP secondary structure sequence.

        """
        return "".join(r.ss for r in self.residues)


class DSSPProtein(object):
    """A single protein parsed from DSSP output (consists of one or more
    chains).

    """

    def __init__(self, pdb_fn):
        """Constructor.

        Args:
            pdb_fn: Path to a PDB file.

        """

        self.dssp = None
        self.chains = {}

        self.parse_dssp(pdb_fn)

        return

    def __getitem__(self, chain_id):
        return self.chains[chain_id]

    def parse_dssp(self, pdb_fn, interchain_bps=False):
        """Parse the output of DSSP when run on a single PDB file.

        Args:
            pdb_fn: Path to a PDB file.
            interchain_bps: True if considering bridge-pairings between chains (optional).

        Returns:
            None

        """

        struct = PDBParser().get_structure(pdb_fn, pdb_fn)
        model = struct[0]
        dssp = DSSP(model, pdb_fn)
        
        # Map each PDB resid to the corresponding DSSP residue information.
        pdb_resid_to_dssp_info = OrderedDict()

        # Map each DSSP index to the corresponding DSSP residue information.
        dssp_index_to_dssp_info = OrderedDict()

        # Begin parsing the DSSP output.
        #
        # TODO:
        # - refactor these DSSP-related code blocks to use BioPython's DSSP
        #   parsing instead.
        with Popen(("dsspcmbi %s" % pdb_fn).split(), stdout=PIPE).stdout as f:

            # Skip the column header line.
            for line in f:
                if line.strip().startswith("#  RESIDUE"):
                    break

            for line in f:
                # Skip chain breaks.
                if line[13] == "!":
                    continue

                pdb_resseq = int(line[5:10])
                pdb_icode = line[10]
                pdb_resid = (pdb_resseq, pdb_icode)
                dssp_index = int(line[:5])
                bp1 = int(line[25:29])
                bp2 = int(line[29:33])

                hb_don_1 = int(line[38:45])
                hb_don_1_energy = float(line[46:50])
                hb_don_1_dssp_index = dssp_index + hb_don_1 if hb_don_1 != 0 else None

                hb_acc_1 = int(line[50:56])
                hb_acc_1_energy = float(line[57:61])
                hb_acc_1_dssp_index = dssp_index + hb_acc_1 if hb_acc_1 != 0 else None

                hb_don_2 = int(line[61:67])
                hb_don_2_energy = float(line[68:72])
                hb_don_2_dssp_index = dssp_index + hb_don_2 if hb_don_2 != 0 else None

                hb_acc_2 = int(line[72:78])
                hb_acc_2_energy = float(line[57:61])
                hb_acc_2_dssp_index = dssp_index + hb_acc_2 if hb_acc_2 != 0 else None

                pdb_resid_to_dssp_info[pdb_resid] = \
                    dssp_index_to_dssp_info[dssp_index] = \
                        { "dssp_index" : dssp_index,
                          "bp1" : bp1, "bp2" : bp2,
                          "hb_don_1_dssp_index" : hb_don_1_dssp_index,
                          "hb_don_1_energy" : hb_don_1_energy,
                          "hb_acc_1_dssp_index" : hb_acc_1_dssp_index,
                          "hb_acc_1_energy" : hb_acc_1_energy,
                          "hb_don_2_dssp_index" : hb_don_2_dssp_index,
                          "hb_don_2_energy" : hb_don_2_energy,
                          "hb_acc_2_dssp_index" : hb_acc_2_dssp_index,
                          "hb_acc_2_energy" : hb_acc_2_energy }

        # Map each DSSP index to its Bio.PDB.Residue object.                
        dssp_index_to_res = {}

        # Update each Bio.PDB.Residue with its DSSP information.
        for val in dssp:
            res = val[0]
            ss = val[1]

            pdb_resid = res.id[1:]
            dssp_info = pdb_resid_to_dssp_info[pdb_resid]

            res.pdb_resid = pdb_resid
            res.dssp_index = dssp_info["dssp_index"]
            res.chain_id = res.get_full_id()[2]
            res.ss = ss

            dssp_index_to_res[res.dssp_index] = res

        raw_chains = defaultdict(list)

        for res in (x[0] for x in dssp):
            pdb_resid = res.id[1:]
            dssp_info = pdb_resid_to_dssp_info[pdb_resid]

            bp1_res = dssp_index_to_res.get(dssp_info["bp1"], None)
            bp2_res = dssp_index_to_res.get(dssp_info["bp2"], None)

            # Ignore interchain bridge partners if required.
            if bp1_res and not interchain_bps:
                if bp1_res.chain_id != res.chain_id:
                    bp1_res = None

            if bp2_res and not interchain_bps:
                if bp2_res.chain_id != res.chain_id:
                    bp2_res = None

            hb_don_1_dssp_index = dssp_info["hb_don_1_dssp_index"]
            hb_don_1_dssp_index = dssp_info["hb_don_1_dssp_index"]

            hb_acc_1_dssp_index = dssp_info["hb_acc_1_dssp_index"]
            hb_acc_1_dssp_index = dssp_info["hb_acc_1_dssp_index"]

            hb_don_2_dssp_index = dssp_info["hb_don_2_dssp_index"]
            hb_don_2_dssp_index = dssp_info["hb_don_2_dssp_index"]

            hb_acc_2_dssp_index = dssp_info["hb_acc_2_dssp_index"]
            hb_acc_2_dssp_index = dssp_info["hb_acc_2_dssp_index"]

            hb_don_1_energy = dssp_info["hb_don_1_energy"]
            hb_don_1_energy = dssp_info["hb_don_1_energy"]

            hb_acc_1_energy = dssp_info["hb_acc_1_energy"]
            hb_acc_1_energy = dssp_info["hb_acc_1_energy"]

            hb_don_2_energy = dssp_info["hb_don_2_energy"]
            hb_don_2_energy = dssp_info["hb_don_2_energy"]

            hb_acc_2_energy = dssp_info["hb_acc_2_energy"]
            hb_acc_2_energy = dssp_info["hb_acc_2_energy"]

            hb_don_1_res = dssp_index_to_res.get(hb_don_1_dssp_index, None)
            hb_don_2_res = dssp_index_to_res.get(hb_don_2_dssp_index, None)

            hb_acc_1_res = dssp_index_to_res.get(hb_acc_1_dssp_index, None)
            hb_acc_2_res = dssp_index_to_res.get(hb_acc_2_dssp_index, None)

            # TODO:
            # - Replace these lines with more elegant code.
            res.bp1 = bp1_res
            res.bp2 = bp2_res

            res.hb_don_1_res = hb_don_1_res
            res.hb_don_2_res = hb_don_2_res

            res.hb_acc_1_res = hb_acc_1_res
            res.hb_acc_2_res = hb_acc_2_res

            res.hb_don_1_dssp_index = hb_don_1_dssp_index
            res.hb_don_2_dssp_index = hb_don_2_dssp_index

            res.hb_acc_1_dssp_index = hb_acc_1_dssp_index
            res.hb_acc_2_dssp_index = hb_acc_2_dssp_index

            res.hb_don_1_energy = hb_don_1_energy
            res.hb_don_2_energy = hb_don_2_energy

            res.hb_acc_1_energy = hb_acc_1_energy
            res.hb_acc_2_energy = hb_acc_2_energy

            raw_chains[res.chain_id].append(res)

        for chain_id, raw_chain in raw_chains.iteritems():
            self.chains[chain_id] = DSSPChain(raw_chain)

        self.dssp = dssp

        return


class Strand(object):
    """A wrapper class around PTNodeStrand representing a single beta-strand.

    """

    def __init__(self, pdb_id, chain_id, row_index, ptnodestrand):
        """Constructor.

        Arguments:
            pdb_id -- unique PDB id of the parent protein (e.g. "1UBQ").
            ptnodestrand -- ptnode.PTNodeStrand object to be wrapped.
            topology_row --

        """

        self.pdb_id = pdb_id
        self.chain_id = chain_id

        self.ptnodestrand = ptnodestrand
        self.align_pos = ptnodestrand.align_pos
        self.reversed = ptnodestrand.reversed
        self.is_barrel_edge = ptnodestrand.barrel_edge
        self.dir_str = "<" if self.reversed else ">"

        self.row_index = row_index

        # the strand_id is a globally-unique strand identifier, global over all
        # other chains, sheets, and strands in the parent protein.
        self.strand_id = "%s%s_%s" % (pdb_id.upper(), ptnodestrand.chainid,
                                      ptnodestrand.nodeid)
        self.strand_num = int(ptnodestrand.nodeid.split("_")[-1])

        # ordinal list of Bio.PDB.Residue objects, use this list for iterating
        # through the strand residues in sequential order..
        self.ord_res_list = []
        self.res_dict = defaultdict(bool)

        self.__build_ordinal_res_list()

        return

    @property
    def ptns(self): return self.ptnodestrand

    @property
    def is_barrel_edge(self): return self.ptns.get_barrel_edge()

    def __len__(self):
        """Get the number of residues in this strand.

        """

        return len(self.ord_res_list)

    def __repr__(self):
        """Get a string representation of the strand (ie. the strand_id).

        """

        return self.strand_id

    def __iter__(self):
        """

        """

        return iter(self.ord_res_list)

    def iter_pdb_resids(self):
        return ( r.pdb_resid for r in self )

    def __hash__(self):
        """Get the hash value of this Strand object as the hash of its
        strand_id.

        """

        return hash(self.strand_id)

    def __build_ordinal_res_list(self):
        """Build the ordinal list of Bio.PDB.Residue objects in this strand.

        Arguments:
            None

        Returns:
            None

        """

        strand_residues = self.ptnodestrand.get_residue_list()

        self.ord_res_list = [ None for _ in strand_residues ]

        for res in strand_residues:
            if self.reversed:
                res.strand_dir = "<"
            else:
                res.strand_dir = ">"

            res.strand_id = self.strand_id
            res.strand_num = self.strand_num
            ord_index = self.ptnodestrand.get_residue_ordinal(res.pdb_resid) - 1
            self.ord_res_list[ord_index] = res
            self.res_dict[res.pdb_resid] = res

        return

    @property
    def direction(self):
        direction = "<" if self.reversed else ">"
        return direction

    def get_node_label(self):
        """
        """

        return self.strand_id

    def is_barrel_edge(self):
        """

        """

        return self.ptnodestrand.get_barrel_edge()

    def is_parallel(self, other):
        return self.ptnodestrand.is_parallel(other.ptnodestrand)


class Sheet(object):
    """A class representing a beta-sheet.

    """

    def __init__(self, pdb_id, chain_id, ptsheet_id, ptnodestrand_grid, ptg,
            chain_res_list):
        """Constructor.

        Arguments:
            

        Returns:
            None

        """

        # sheet meta-data.
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.sheet_id = self.make_sheet_id(pdb_id, chain_id, ptsheet_id)

        # set to True if this sheet is barreled, otherwise set to False.
        #
        # NOTE: if this sheet is barrelled, then ptgraph2 would have set the
        # PTNodeStrand.barrel_edge attribute to two strands in this sheet to
        # True.
        self.is_barrel = False

        # set to True if this sheet contains bifurcations i.e. two strands on
        # the same side of a reference strand.
        self.is_bifurcated = False

        # keep a list of all the barreled strands (Strand objects) in this
        # sheet. A barrelled strand is one of the two strands (staves of the
        # barrel) that have been "opened" up to create a planar representation
        # of the beta-sheet. The choice of which strands to open up was
        # arbitrarily made by ptgraph2.
        self.barreled_strands = []

        # strand pairing graph where nodes are labelled according to the
        # self.get_node_label() function.
        # edges represent bridge-pairings between strands.
        # NOTE: interchain bridge pairings are ignored.
        self.strand_G = None

        # residue pairing graph (see documentation in self.__build_graphs(...)
        # for a more detailed explanation).
        self.res_G = None

        # save a reference to the ptnodestrand_grid object that this sheet is
        # made from.

        # a two-nested list containing the rank and order of each strand in a
        # beta-sheet plane.
        self.topology = []

        # lookup Strand objects by their PTNode.nodeid id (NOTE: this id isn't
        # globally unique).
        self.ptnode_id_to_strand = {}

        # lookup Strand objects by their strand_id's (Strand.strand_id).
        self.strand_dict = defaultdict(bool)

        # lookup Strand objects by their ordinal numbers
        # (Strand.ord_strand_num).
        self.strand_ord_dict = OrderedDict()

        self.max_strand_length = 0

        self.res_dict = {}

        # list of Strand objects in order of their appearance in their parent
        # chain (ordinal).
        self.strand_list = []

        self.chain_res_by_resid = OrderedDict()
        self.chain_res_by_ord = OrderedDict()
        self.numeric_topology = []

        for res in chain_res_list:
            self.chain_res_by_resid[ res.pdb_resid ] = res
            self.chain_res_by_ord[ res.ord_index ] = res

        self.ptg = ptg
        self.__build_dicts_and_topology(ptnodestrand_grid)
        self.__build_graphs(ptg)

        return

    def iter_chain_res(self):
        return self.chain_res_by_ord.itervalues()

    def __len__(self):
        """Get the number of strands in this sheet.

        """

        return len(self.strand_dict)

    def __repr__(self):
        """Get the sheet_id as the string representation of this sheet.

        """

        return self.sheet_id

    def __iter__(self):
        """Get an iterator over the strands in the sheet.

        """

        return iter(self.strand_list)

    def __build_dicts_and_topology(self, ptnodestrand_grid):
        """Build the dictionaries and native strand topology associated with
        this Sheet.

        Arguments:
            ptnodestrand_grid --

        Returns:
            None

        """

        # build the strand grid and strand dict.
        for row_index, row in enumerate(ptnodestrand_grid):
            topology_row = []
            numeric_topology_row = []
            for ptnodestrand in row:
                # perform a sanity check for strands with different chain id's
                # to this sheet.
                assert ptnodestrand.chainid == self.chain_id, \
                        ("Attempting to insert a strand from a different chain"
                         "this sheet.")

                strand = Strand(self.pdb_id, self.chain_id, row_index,
                                ptnodestrand)

                if len(strand) > self.max_strand_length:
                    self.max_strand_length = len(strand)

                self.strand_ord_dict[strand.strand_num] = strand
                self.strand_dict[strand.strand_id] = strand
                self.strand_list.append(strand)
                self.ptnode_id_to_strand[ptnodestrand.nodeid] = strand

                if ptnodestrand.barrel_edge:
                    self.is_barrel = True
                    self.barreled_strands.append(strand)

                topology_row.append(strand)
                numeric_topology_row.append(strand.strand_num)

                # assign the native topology row index to each Residue object.
                for res in strand:
                    res.direction = strand.dir_str
                    res.row_index = row_index
                    res.ord_index = \
                        self.ptg.pdb_resid_dict[(res.chain_id, res.pdb_resid)]
                    self.res_dict[res.pdb_resid] = res

            # sort the topology row by strand-offsets (Strand.get_align_pos()).
            topology_row.sort(key=lambda strand : -strand.align_pos)

            self.topology.append(topology_row)
            self.numeric_topology.append(numeric_topology_row)

        return

    def __build_graphs(self, ptg):
        """Build the strand-pairing and residue graphs.

        <self.res_G> is the residue-pairing graph illustrating the inter-residue
        relationships within the protein chain/sheet. These inter-residue
        relationships may include, but is not limited to, peptide bonds (implied
        by the amino acid sequence order, COVALENT interaction), and
        non-covalent interactions such as bridge-pairings
        between spacially-adjacent residues within beta-sheets, disulphide bonds
        between spacially-close cysteine-cysteine pairs, and hydrogen-bonds
        (stabilising the beta-sheet structure).
        
        In the future we want to add extra relationships such as salt-bridges
        (see http://bit.ly/dDv9nW, wikipedia) between oppositely-charged
        side-chains that are also hydrogen-bonded e.g. glutamic acid (-) ::
        lysine (+).

        Returns:
            None

        Raises:
            SelfBridgeException if the sheet contains a bridge between residues
            in the same strand.

        """

        self.strand_G = networkx.Graph()
        self.res_G = networkx.MultiGraph()

        # add the nodes to the res_G.
        for res in self.res_dict.itervalues():
            self.res_G.add_node(res.pdb_resid,
                                label=res.single_aa,
                                res=res)

        # add bridge edges between strand vertices.
        for src_strand in self:
            # add peptide edges between sequential residues in each strand.
            for res_i in xrange(len(src_strand) - 1):
                src_res = src_strand.ord_res_list[res_i]
                des_res = src_strand.ord_res_list[res_i + 1]

                # ensure that edges are only being added between nodes belonging
                # to this sheet.
                assert(self.res_G.has_node(src_res.pdb_resid))
                assert(self.res_G.has_node(des_res.pdb_resid))

                self.res_G.add_edge(src_res.pdb_resid,
                                    des_res.pdb_resid,
                                    label="peptide")

            ptbridge_list = src_strand.ptnodestrand.get_bridge_list()

            for des_ptnode, direction, side in ptbridge_list:
                des_strand = self.ptnode_id_to_strand[des_ptnode.nodeid]
                self.strand_G.add_edge(src_strand.get_node_label(),
                                       des_strand.get_node_label(),
                                       side=side, label="bridge")

        # determine if this sheet is bifurcated
        for rank in self.topology:
            if len(rank) > 1:
                self.is_bifurcated = True
                break

        bridge_edges = []

        # go through all the residue bridges in ptg.secstruct.bridgeres_list and
        # only add bridge edges between residues in this sheet.
        for src_chain_id, src_pdb_resseq, des_chain_id, des_pdb_resseq, bdir in\
            ptg.secstruct.bridgeres_list:

            # ignore inter-chain bridges.
            if not (src_chain_id == des_chain_id == self.chain_id):
                continue

            # ignore inter-sheet bridges.
            if src_pdb_resseq not in self.res_dict or \
               des_pdb_resseq not in self.res_dict:
                continue

            src_res = self.res_dict[src_pdb_resseq]
            des_res = self.res_dict[des_pdb_resseq]

            # ensure that edges are only being added between nodes belonging
            # to this sheet.
            assert(self.res_G.has_node(src_res.pdb_resid))
            assert(self.res_G.has_node(des_res.pdb_resid))

            # if there is a bridge edge between a pair of residues within the
            # same strand then dssp/stride/pdb_secstruct secondary assignment
            # methods have missed a beta-turn classification ('T') thus not
            # delineating the strand pairs from each other. Due to the inherent
            # inaccuracies of said secondary assignment methods, this
            # pathological case can conceivably occur alot.
            if src_res.strand_id == des_res.strand_id:
                raise SelfBridgeException, \
                        ("Bridge between residues %s and %s"
                         " occur within the same strand (%s)" %
                         (src_res.pdb_resid, des_res.pdb_resid,
                          src_res.strand_id))

            edge = (src_res.pdb_resid, des_res.pdb_resid)

            if edge in bridge_edges or tuple(reversed(edge)) in bridge_edges:
                continue

            bridge_edges.append(edge)

            self.res_G.add_edge(src_res.pdb_resid,
                                des_res.pdb_resid,
                                label="bridge", bdir=bdir)

        assert self.res_G.number_of_nodes() > 0, \
                ("no residue nodes discovered in %s" % self.sheet_id)

        return

    @classmethod
    def make_sheet_id(cls, pdb_id, chain_id, ptsheet_id):
        """

        """

        return "sheet-%s-%s-%s" % (pdb_id.lower(), chain_id, ptsheet_id)

    def build_sheet_matrix(self):
        """Build a two-dimensional matrix of residue node id's as the planar
        representation of this beta-sheet.

        This function uses the griddify.py script written by Graeme Gange
        <ggange@csse.unimelb.edu.au> that generates this matrix from a list of
        nodes and bridge edges.

        Returns:
            a SheetMatrix object representation of this Sheet object.

        """

        try:
            sheet_matrix = SheetMatrix(self)
        except:
            raise BadSheetMatrixException

        return sheet_matrix

    def generate_string(self):
        """Generate the one-dimensional sequence representation of this sheet
        matrix.

        

        """

        sheet_matrix = self.build_sheet_matrix()

        for col in len(sheet_matrix[0]):
            for row in xrange(len(sheet_matrix)):
                if sheet_matrix[row][col] == None:
                    # pseudo node
                    pass
                else:
                    # real node.
                    pass

        return

    def iter_residues(self):
        """
        """

        for strand in self:
            for res in strand:
                yield res
        return


class Chain(object):
    """A class representing a protein chain.

    """

    def __init__(self, pdb_id, chain_id, res_list, ptg):
        """Constructor.

        Arguments:
            pdb_id --
            chain_id --
            res_list -- list of Residue objects ordered by their ordinal
            indices in the Chain.

        """

        self.pdb_id = pdb_id
        self.chain_id = chain_id
        
        ## These contain Bio.PDB.Residue.Residue objects.
        self.res_list = res_list
        self.res_dict = {} # Lookup residues by their pdb_resid's.

        ## This contains betapy.Sheet objects
        self.sheets_dict = {} # lookup sheets by their sheet ids.

        self.__build_res_dict()
        self.__build_secondary_structures(ptg)
        
        return

    def __len__(self):
        return len(self.res_list)

    def __iter__(self):
        """Return an iterator over the residues in this Chain.

        Returns:
            an iterator over the residues in this Chain.

        """

        return iter(self.res_list)

    def __build_res_dict(self):
        """Build the residue lookup dictionary.

        Arguments:
            None

        Returns:
            None

        """

        for res in self.res_list:
            self.res_dict[res.pdb_resid] = res

        return

    def __build_secondary_structures(self, ptg):
        """

        """

        # iterator through all the sheets in the ptgraph object.
        ptsheet_iter = ptg.sheet_strandlists_dict.iteritems()

        for ptsheet_id, ptnodestrand_grid in ptsheet_iter:
            not_this_chain = False

            for row in ptnodestrand_grid:
                for ptnodestrand in row:
                    if ptnodestrand.chainid != self.chain_id:
                        not_this_chain = True
                        break
                if not_this_chain:
                    break

            # ignore a ptnodestrand_grid if it contains strands not of this
            # chain (e.g. inter-chain sheets).
            if not_this_chain:
                continue

            sheet_id = Sheet.make_sheet_id(self.pdb_id, self.chain_id,
                                           ptsheet_id)

            try:
                sheet = Sheet(self.pdb_id, self.chain_id, ptsheet_id,
                              ptnodestrand_grid, ptg, self.res_list)
            except SelfBridgeException:
                # sheet contains self bridges, do nothing, don't add the sheet
                # to this chain.
                sys.stderr.write("WARNING: IGNORING SHEET ID (SELF-BRIDGE): %s\n" % sheet_id)
            except:
                sys.stderr.write("WARNING: IGNORING SHEET ID: %s\n" % sheet_id)
            else:
                self.sheets_dict[sheet.sheet_id] = sheet

        return

    def sheets_iter(self):
        """Get an iterator over the Sheet objects of this Chain.

        """

        return self.sheets_dict.itervalues()


class AbstractMatrix(list):
    """An abstract class representing an N x M beta-sheet where N is the number
    of rows and M is the number of columns. Each row may contain more than one
    beta-strand.

    """

    def __init__(self, data):
        """Constructor

        """
        list.__init__(self, data)

        return

    def __find_peptide_neighbour(self, coord, direction):
        """Find the col-index of the next residue that is peptide bonded to the
        one at self[row][col] in a given direction.

        """

        row, col = coord
        new_col = col

        if direction == "left":
            delta = -1
        else:
            delta = 1

        while self[row][new_col] == '-':
            new_col += delta

        return new_col

    def is_peptide_bond(self, row, col):
        return self[row][col] == '-'

    def get_single_aa(self, i, j):
        """Get the single-letter amino acid representation of the entry at
        self[i][j].

        """

        raise NotImplementedError("This function has not been implemented.")

        return

    def get_csv_string(self):
        """Get a csv string representation of this SheetMatrix.

        """

        if self.is_barrel:
            barrel_str = 't'
        else:
            barrel_str = 'f'

        if self.is_bifurcated:
            bifurcated_str = 't'
        else:
            bifurcated_str = 'f'

        is_bulged = False

        mat_rows = []
        max_strand_length = 0

        new_topology = []

        strand_directions = OrderedDict()

        for group in self._sheet.topology:
            new_group = []
            for strand in group:
                if strand.reversed:
                    direction_str = "<"
                else:
                    direction_str = ">"
                strand.dir_str = direction_str
                new_group.append(str(strand.strand_num))
            new_topology.append(":".join(new_group))

        strands_dict = OrderedDict()

        for strand_num, strand in self._sheet.strand_ord_dict.iteritems():
            strand_header = "%d %s" % (strand_num, strand.dir_str)
            strand_str = strand_header + " " + " ".join("{r.single_aa}:{r.pdb_resid}:{r.ord_index}".format(r=r) for r in strand)
            strands_dict[strand.strand_num] = strand_str

        for i, row in enumerate(self):
            row_list = []
            for j, val in enumerate(row):
                if val == None:
                    row_list.append('.')
                elif val == '-':
                    row_list.append('-')
                    is_bulged = True
                else:
                    _, res, _ = self.get_single_aa(i, j)
                    row_list.append(res.get_csv_str())

            strand_str = " ".join(row_list)
            mat_rows.append(strand_str)

        row_vals = [ self.sheet_id, self.pdb_id, self.chain_id ]
        row_vals.append(" ".join("%s" % v for v in new_topology))
        row_vals.append(";".join("%s" % v for v in strands_dict.itervalues()))
        row_vals.append(";".join(mat_rows))
        # row_vals.append(" ".join(r.get_csv_str() for r in iter_chain_res))

        return ";".join( row_vals )


class SheetMatrix(AbstractMatrix):
    """This class also wraps a Sheet object.

    """

    def __init__(self, sheet):
        """Constructor.

        Arguments:
            data -- a "raw" two-dimensional list representation of the sheet
            matrix, this should be generated by griddify.build_grid(...).

        Returns:
            None

        """

        self._sheet = sheet
        self.chain_id = sheet.chain_id.upper()
        self.topology = [ [ x.strand_num for x in s ] for s in sheet.topology ]
        self.orients = {}
        self.strand_num_dict = sheet.strand_ord_dict
        self.strand_nums = self.strand_num_dict.keys()

        for src, des in sheet.strand_G.edges_iter():
            src_strand = self.strand_dict[src]
            des_strand = self.strand_dict[des]

            for res in src_strand:
                assert(res.strand_num == src_strand.strand_num)

            for res in des_strand:
                assert(res.strand_num == des_strand.strand_num)

            key = frozenset((src_strand.strand_num, des_strand.strand_num))
            val = src_strand.is_parallel(des_strand)
            self.orients[key] = val

        self.is_bulged = False
        self.__build_sheet_matrix(sheet)

        return

    def __repr2__(self):
        """Get a string representation of this class

        """

        if self.is_barrel:
            barrel_str = 't'
        else:
            barrel_str = 'f'

        if self.is_bifurcated:
            bifurcated_str = 't'
        else:
            bifurcated_str = 'f'

        n_edges = 0

        # count the number of edges.
        for i in xrange(len(self)):
            for j in xrange(len(self[i])):
                if self[i][j] == None or self[i][j] == '-':
                    continue

                if i < len(self) - 1:
                    # check down neighbour.
                    if self[i+1][j] == None or self[i+1][j] == '-':
                        pass
                    elif self.get_single_aa(i + 1, j)[0] in aa1:
                        n_edges += 1
                else:
                    # don't check down neighbour.
                    pass

                if j < len(self[i]) - 1:
                    # check right neighbour.
                    next_col = j + 1

                    while self[i][next_col] == '-':
                        next_col += 1

                    if self[i][next_col] == None:
                        pass
                    elif self.get_single_aa(i, next_col)[0] in aa1:
                        n_edges += 1
                else:
                    # don't check right neighbour.
                    pass
        
        n_nodes = self.res_G.number_of_nodes()

        if self.is_barrel:
            n_nodes *= 2

        n_strands = len(self.strand_nums)
        buf = "%s;%d;%d;%d;%s;%s;%d;" % (self.sheet_id, len(self), len(self[0]),
                                         n_strands, bifurcated_str, barrel_str, n_nodes)

        for i, row in enumerate(self):
            row_list = []
            for j, val in enumerate(row):
                if val == None:
                    row_list.append('.')
                elif val == '-':
                    row_list.append('-')
                else:
                    row_list.append(self.get_single_aa(i, j)[0])

            buf += ''.join(row_list) + ','

        return buf[:-1]

    def iter_strand_pairs(self):
        return ( (int(x.split("_")[-1]),
                  int(y.split("_")[-1]))
                 for (x, y) in self._sheet.strand_G.edges_iter() )

    def __getattr__(self, attr):
        """Redirect attribute references to self._sheet if not found in
        SheetMatrix.
        
        Arguments:
            attr -- name of the attribute.

        Returns:
            a reference to the requested attribute.

        """

        return getattr(self._sheet, attr)
    



    @property
    def sheet(self):
        return self._sheet

    @property
    def is_bifurcated(self):
        return self.sheet.is_bifurcated

    @property
    def is_barrel(self):
        return self.sheet.is_barrel

    @property
    def n_strands(self):
        return len(self.strand_num_dict)

    @property
    def max_strand_length(self):
        return self.sheet.max_strand_length




    def __add_peptide_nodes(self, raw_matrix):
        """

        """

        for row in raw_matrix:
            for i in xrange(len(row) - 1):
                if row[i] == None or row[i + 1] != None:
                    continue

                j = i + 1

                # find the row neighbour residue in the row, skipping over dummy
                # nodes.
                while row[j] == None and j < len(row) - 1:
                    j += 1

                # search the next row if there are no neighbours on this row.
                if row[j] == None:
                    continue

                # check if row[i] and row[j] are on the same strand.
                i_strand_id = self.get_res_by_nodeid(row[i]).strand_id
                j_strand_id = self.get_res_by_nodeid(row[j]).strand_id

                if i_strand_id == j_strand_id:
                    # place "stretched" peptide bond entries between row[i] and
                    # row[j]:
                    for k in xrange(i + 1, j):
                        row[k] = '-'

        return

    def __add_residue_objects(self, raw_matrix):
        """

        """

        # replace node ids with Residue objects in the raw sheet matrix.
        for row in raw_matrix:
            for i in xrange(len(row)):
                if row[i] == None:
                    continue
                elif row[i] == '-':
                    self.is_bulged = True
                    continue
                else:
                    row[i] = self.get_res_by_nodeid(row[i])

        return

    def __are_bridges_satisfied(self, raw_matrix, bridges):
        """

        """

        # remove the bridge directions (e.g. 'N'/'P') from the bridges.
        unsatisfied_pairs = [tuple(sorted(p)) for ((p), _) in bridges]

        # Check if all bridge edges have been satisfied by the raw_matrix.
        for row in xrange(len(raw_matrix) - 1):
            for col in xrange(len(raw_matrix[row])):
                if raw_matrix[row][col] == None:
                    continue

                curr_pair = tuple(sorted((raw_matrix[row][col], 
                                          raw_matrix[row + 1][col])))

                if curr_pair in unsatisfied_pairs:
                    unsatisfied_pairs.remove(curr_pair)

        if len(unsatisfied_pairs) != 0:
            raise BridgesNotSatisfiedException

        return

    def __build_sheet_matrix(self, sheet):
        """

        """

        raw_matrix = []
        all_res_nodes = []

        # Construct a three-nested list of residue nodes from its strands.
        # converts this:
        #
        #   sheet.topology == [[<Strand>, <Strand>], [<Strand>, ...]]
        #
        # Into this:
        #
        #   raw_matrix == [[[1, 2, 3], [10, 11, 12]], [[20, 19, 18], ...]]
        #
        for i, row in enumerate(sheet.topology):
            mat_row = []
            for strand in row:
                res_nodes = [res.pdb_resid for res in strand.ord_res_list]
                if strand.reversed:
                    res_nodes.reverse()
                all_res_nodes.extend(res_nodes)
                mat_row.append(res_nodes)
            raw_matrix.append(mat_row)

        if self.is_barrel:
            dup_raw_matrix = copy.deepcopy(raw_matrix)

            for row in dup_raw_matrix:
                for i in xrange(len(row)):
                    for j in xrange(len(row[i])):
                        row[i][j] = '^' + row[i][j]
                        all_res_nodes.append(row[i][j])

            raw_matrix.extend(dup_raw_matrix)

        # Find all the bridges between residues in this sheet.
        bridges = []

        for src, des, data in sheet.res_G.edges_iter(data=True):
            if data["label"] != "bridge":
                continue

            src_strand = \
                sheet.strand_dict[sheet.res_G.node[src]["res"].strand_id]
            des_strand = \
                sheet.strand_dict[sheet.res_G.node[des]["res"].strand_id]

            if self.is_barrel:
                # Exclude bridges between residues in strands that were chosen
                # as a barrel "flattening point".
                if src_strand.is_barrel_edge and des_strand.is_barrel_edge:
                    if src_strand.row_index < des_strand.row_index:
                        src = '^' + src
                    else:
                        des = '^' + des

                    bridges.append(((src, des), data["bdir"]))
                    continue

                # Store the duplicate bridges.
                bridges.append((('^' + src, '^' + des), data["bdir"]))

            bridges.append(((src, des), data["bdir"]))

        if len(bridges) == 0:
            return None

        raw_sheet_matrix = build_grid(all_res_nodes, raw_matrix, bridges)

        self.__are_bridges_satisfied(raw_sheet_matrix, bridges)
        self.__add_peptide_nodes(raw_sheet_matrix)
        self.__add_residue_objects(raw_sheet_matrix)

        AbstractMatrix.__init__(self, raw_sheet_matrix)

        return

    def get_res_by_nodeid(self, node_id):
        """Get a Residue object inside the SheetMatrix by its node_id.

        """

        if node_id[0] == '^':
            node_id = node_id[1:]

        return self.res_dict[node_id]

    def get_single_aa(self, i, j):
        """Get the single-letter amino acid representation of the entry at
        self[i][j].

        """

        curr_strand = self.strand_dict[self[i][j].strand_id]
        return self[i][j].single_aa, self[i][j], curr_strand

    def get_resids_repr(self, basic=False, halfbarrels=False):
        """Get a string representation of this class

        """

        if self.is_barrel:
            barrel_str = 't'
        else:
            barrel_str = 'f'

        n_edges = 0

        # count the number of edges.
        for i in xrange(len(self)):
            for j in xrange(len(self[i])):
                if self[i][j] == None or self[i][j] == '-':
                    continue

                if i < len(self) - 1:
                    # check down neighbour.
                    if self[i+1][j] == None or self[i+1][j] == '-':
                        pass
                    else:
                        aa, bpy_res, _ = self.get_single_aa(i+1,j)
                        if aa in aa1:
                            n_edges += 1
                else:
                    # don't check down neighbour.
                    pass

                if j < len(self[i]) - 1:
                    # check right neighbour.
                    next_col = j + 1

                    while self[i][next_col] == '-':
                        next_col += 1

                    if self[i][next_col] == None:
                        pass
                    else:
                        aa, bpy_res, _ = self.get_single_aa(i, next_col)
                        if aa in aa1:
                            n_edges += 1
                else:
                    # don't check right neighbour.
                    pass
        
        n_nodes = self.res_G.number_of_nodes()

        if self.is_barrel:
            n_nodes *= 2

        if basic:
            buf = "%s:" % self.sheet_id
        else:
            buf = "%s:0:%d:%d:%s:%d:" % \
                    (self.sheet_id, len(self), len(self[0]),
                     barrel_str, n_nodes)

        rows = enumerate(self)

        if halfbarrels:
            if self.is_barrel:
                rows = ( (i,r) for i,r in enumerate(self) if i < (len(self) / 2) )

        for i, row in rows:
            row_list = []
            for j, val in enumerate(row):
                if val == None:
                    row_list.append('.')
                elif val == '-':
                    row_list.append('-')
                else:
                    aa, _, _ = self.get_single_aa(i,j)
                    row_list.append(aa + "#" + self[i][j].pdb_resid + "#")

            buf += '|'.join(row_list) + ','

        return buf[:-1]

    def get_stats_parseable_str(self):
        """
        """
        barrel_str = 't' if self.is_barrel else 'f'

        G = self.strand_G

        for n in G.nodes_iter():
            if G.degree(n) > 2:
                bifurcated_str = 't'
            else:
                bifurcated_str = 'f'

        buf_fields = [self.sheet_id, barrel_str, bifurcated_str]

        # get the strand orient pairs.
        orient_list = []

        for (src, des), is_parallel in self.orients.iteritems():
            orient = 'p' if is_parallel else 'n'  
            orient_list.append("%d^%d^%s" % (src, des, orient))

        buf_fields.append("+".join(orient_list))

        for i, row in enumerate(self):
            row_list = []
            for j, val in enumerate(row):
                if val == None:
                    row_list.append('.')
                elif val == '-':
                    row_list.append('-')
                else:
                    res = self[i][j]
                    res_str = "%s^%s" % (res.pdb_resid, res.strand_num)
                    row_list.append(res_str)

            buf_fields.append("|".join(row_list))

        return ":".join(buf_fields)
    
    def get_bmat_string(self, pdbpath):
        """Get a string representation of this class

        """

        if self.is_barrel:
            barrel_str = 't'
        else:
            barrel_str = 'f'

        if self.is_bifurcated:
            bifurcated_str = 't'
        else:
            bifurcated_str = 'f'

        n_edges = 0

        # count the number of edges.
        for i in xrange(len(self)):
            for j in xrange(len(self[i])):
                if self[i][j] == None or self[i][j] == '-':
                    continue

                if i < len(self) - 1:
                    # check down neighbour.
                    if self[i+1][j] == None or self[i+1][j] == '-':
                        pass
                    elif self.get_single_aa(i + 1, j)[0] in aa1:
                        n_edges += 1
                else:
                    # don't check down neighbour.
                    pass

                if j < len(self[i]) - 1:
                    # check right neighbour.
                    next_col = j + 1

                    while self[i][next_col] == '-':
                        next_col += 1

                    if self[i][next_col] == None:
                        pass
                    elif self.get_single_aa(i, next_col)[0] in aa1:
                        n_edges += 1
                else:
                    # don't check right neighbour.
                    pass
        
        n_nodes = self.res_G.number_of_nodes()

        if self.is_barrel:
            n_nodes *= 2

        #
        # Parse the PDB file header in order to get the chain and sequence
        # metadata.
        #

        pdb_id = self.sheet_id.split("-")[1].lower()
        chain_id = self.chain_id.upper()
        pdb_fn = os.path.join(os.path.expanduser(pdbpath),
                pdb_id[1:3], "pdb%s.ent.gz" % pdb_id)

        parser = PDBParser()

        with gzip.open(pdb_fn) as f: 
            structure = parser.get_structure(self.pdb_id, f)

        header = parser.get_header()
        mol_id = None
        molecule_name = "N/A"
        organism_common_name = "N/A"
        organism_scientific_name = "N/A"
        compound_dict = header.get("compound", {})
        source_dict = header.get("source", {})

        for _mol_id, mol_dict in compound_dict.iteritems():
            chains = mol_dict.get("chain", "")
            if chain_id.lower() in chains or chain_id.upper() in chains:
                mol_id = _mol_id
                molecule_name = mol_dict.get("molecule", "N/A")
                break

        if source_dict:
            if mol_id in source_dict:
                organism_common_name = source_dict[mol_id].get("organism_common", "N/A")
                organism_scientific_name = source_dict[mol_id].get("organism_scientific", "N/A")

        n_strands = len(self.strand_nums)
        buf = "%s^%d^%d^%d^%s^%s^%d^%s^%s^%s^" % \
                (self.sheet_id, len(self), len(self[0]),
                 n_strands, bifurcated_str, barrel_str, n_nodes,
                 molecule_name, organism_common_name, organism_scientific_name)

        for i, row in enumerate(self):
            row_list = []
            for j, val in enumerate(row):
                if val == None:
                    row_list.append('.')
                elif val == '-':
                    row_list.append('-')
                else:
                    row_list.append(self.get_single_aa(i, j)[0])

            buf += ''.join(row_list) + ','

        return buf[:-1]

    def iter_strand_pair_json_strings(self):
        """

        """

        for pair in self.iter_strand_pairs():

            print ">>", pair

            yield self.get_json_string(pair)

        return

    def iter_bps(self, duplicates=False, offset=0):
        """

        Returns:
            iterator, where each value is of the form:
                <res_i AA>
                <res_j AA>
                <orientation i.e. "p" | "ap">
                <res_i X1 angle>
                <res_j X1 angle>
                <offset e.g. i, j + offset> 
                <res_i matrix coord>
                <res_j matrix coord>

        """

        NON_RESIDUES = set([None, "-"])

        def get_orient(res_a, res_b):
            if res_a.direction == res_b.direction:
                orient = "p"
            else:
                orient = "ap"

            return orient

        nrows = len(self)
        ncols = len(self[0])

        for i in xrange(nrows):
            if (not duplicates) and (i == nrows - 1):
                break

            for j in xrange(ncols):
                if self[i][j] == None or self[i][j] == '-':
                    continue

                for off in xrange(-offset, offset + 1):
                    if not (0 <= (j + off) < ncols):
                        continue

                    res_a = self[i][j]

                    if i > 0 and duplicates:
                        if self[i-1][j+off] not in NON_RESIDUES:
                            res_b = self[i-1][j+off]
                            orient = get_orient(res_a, res_b)

                            yield (res_a.single_aa,
                                   res_b.single_aa,
                                   res_a.X1_angle,
                                   res_b.X1_angle,
                                   orient,
                                   off,
                                   (i, j),
                                   (i - 1, j + off))

                    if i < nrows - 1 and (self[i+1][j+off] not in NON_RESIDUES):
                        res_b = self[i+1][j+off]
                        orient = get_orient(res_a, res_b)

                        yield (res_a.single_aa,
                               res_b.single_aa,
                               res_a.X1_angle,
                               res_b.X1_angle,
                               orient,
                               off,
                               (i, j),
                               (i + 1, j + off))

        return

    def iter_residues(self):
        for res in self.sheet.iter_residues():
            yield res

    def iter_strands(self):
        """

        """

        for strand in self.sheet:
            yield strand

        return


class Protein(DSSPProtein):
    """
    """

    def __init__(self, pdb_fn):
        """Constructor.

        """

        dssp.DSSPProtein.__init__(self, pdb_fn)

        self.pdb_fn = pdb_fn
        self.ptg = make_graphs(pdb_fn, "none", "dssp")[0]

        self.ptg.build_constraints()
        self.pdb_id = self.ptg.pdb_struct.id.upper()

        self.sheets_by_chainid = defaultdict(list)
        self.sheets_by_sheetid = OrderedDict()
        self.beta_matrices = OrderedDict()

        self._build_sheets()
        self._build_beta_matrices()

        return

    @property
    def sheets(self):
        return self.sheets_by_sheetid.values()

    def _build_sheets(self):
        """
        """

        def make_sheet_id(ptsheet_id):
            return "sheet-%s-%s-%s" % (self.pdb_id.lower(), chain_id, ptsheet_id)

        for chain_id, chain in self.chains.iteritems():
            ptsheet_iter = self.ptg.sheet_strandlists_dict.iteritems()

            for ptsheet_id, ptnodestrand_grid in ptsheet_iter:
                not_this_chain = False

                for row in ptnodestrand_grid:
                    for ptnodestrand in row:
                        if ptnodestrand.chainid != chain_id:
                            not_this_chain = True
                            break

                # ignore a ptnodestrand_grid if it contains strands not of this
                # chain (e.g. inter-chain sheets).
                if not_this_chain:
                    continue

                sheet_id = "sheet-%s-%s-%s" % (self.pdb_id.lower(), chain_id,
                                               ptsheet_id)

                sheet = Sheet(self.pdb_id, chain_id, ptsheet_id,
                              ptnodestrand_grid, self.ptg, chain)

                self.sheets_by_sheetid[sheet.sheet_id] = sheet
                self.sheets_by_chainid[chain_id].append(sheet)

        return

    def _build_beta_matrices(self):
        """
        """

        # Wrap each SheetMatrix using a BetaMatrix.
        for sheet in self.iter_sheets():
            try:
                chain = self[sheet.chain_id]
                beta_matrix = BetaMatrix(sheet.build_sheet_matrix(), chain)
            except:
                sys.stderr.write("WARNING: Could not build BetaMatrix from %s, skipping...\n" % sheet.sheet_id)
                sys.stderr.write(os.linesep)
            else:
                if beta_matrix.n_strands < 2:
                    sys.stderr.write("WARNING: %s contains fewer than 2 strands, skipping...\n" % beta_matrix.sheet_id)
                    continue

                self.beta_matrices[sheet.sheet_id] = beta_matrix

        return

    def __repr__(self):
        return self.pdb_id

    def iter_sheets(self):
        return self.sheets_by_sheetid.itervalues()

    def iter_beta_matrices(self, chain_id=None):
        """

        Arguments:
            chain_id -- chain identifier (default=None).

        """

        for beta_matrix in self.beta_matrices.itervalues():
            if chain_id:
                if beta_matrix.chain_id == chain_id:
                    yield beta_matrix
            else:
                yield beta_matrix

        return


class BetaMatrix(object):
    """This class wraps the betapy.SheetMatrix class.

    """

    def __init__(self, sheet_matrix, chain):
        """Constructor.

        Arguments:
            sheet_matrix -- betapy.SheetMatrix object.
            chain -- dssp.Chain object.

        """

        self._sheet_matrix = sheet_matrix
        self._init_residues(chain)

        return

    @property
    def pdb_id(self):
        return self.sheet.pdb_id

    @property
    def sheet_id(self):
        return self.sheet_matrix.sheet_id

    @property
    def chain_id(self):
        return self.sheet_matrix.chain_id

    @property
    def aa_mat(self):
        """
        """

        def f_single_aa(res):
            if res == None:
                return '.'
            elif type(res) == Bio.PDB.Residue.Residue:
                return res.single_aa
            else:
                return res

        f = numpy.vectorize(f_single_aa)
        mat = f(numpy.array(self.sheet_matrix))

        return mat

    @property
    def pdb_resid_mat(self):
        """
        """

        def f_pdb_resid(res):
            if res == None:
                return '.'
            elif type(res) == Bio.PDB.Residue.Residue:
                return res.pdb_resid
            else:
                return res

        f = numpy.vectorize(f_pdb_resid)
        mat = f(numpy.array(self.sheet_matrix))

        return mat

    @property
    def sheet_matrix(self): return self._sheet_matrix

    @property
    def topology(self): return self.sheet_matrix.topology

    @property
    def is_bifurcated(self): return self.sheet_matrix.is_bifurcated

    @property
    def is_barrel(self): return self.sheet_matrix.is_barrel

    @property
    def is_barreled(self): return self.is_barrel

    @property
    def is_bulged(self): return self.sheet_matrix.is_bulged

    @property
    def n_strands(self): return self.sheet_matrix.n_strands

    @property
    def max_strand_length(self): return self.sheet_matrix.max_strand_length

    @property
    def sheet(self): return self.sheet_matrix.sheet
    
    @property
    def strands(self): return self.sheet.strand_ord_dict

    @property
    def ptg(self): return self.sheet.ptg

    def _init_residues(self, chain):
        """Copy the bridge- and hydrogen-bond partner information from
        the equivalent residues of the DSSPChain object.

        """

        for res in self.iter_residues():
            dssp_res = chain[res.pdb_resid]

            res.dssp_index = dssp_res.dssp_index
            res.ss = dssp_res.ss

            res.acc = dssp_res.acc
            res.phi = dssp_res.phi
            res.psi = dssp_res.psi

            res.bp1 = dssp_res.bp1
            res.bp2 = dssp_res.bp2

            res.hb_don_1_res = dssp_res.hb_don_1_res
            res.hb_don_2_res = dssp_res.hb_don_2_res

            res.hb_acc_1_res = dssp_res.hb_acc_1_res
            res.hb_acc_2_res = dssp_res.hb_acc_2_res

            res.hb_don_1_energy = dssp_res.hb_don_1_energy
            res.hb_don_2_energy = dssp_res.hb_don_2_energy

            res.hb_acc_1_energy = dssp_res.hb_acc_1_energy
            res.hb_acc_2_energy = dssp_res.hb_acc_2_energy

            res.hb_don_1_dssp_index = dssp_res.hb_don_1_dssp_index
            res.hb_don_2_dssp_index = dssp_res.hb_don_2_dssp_index

            res.hb_acc_1_dssp_index = dssp_res.hb_acc_1_dssp_index
            res.hb_acc_2_dssp_index = dssp_res.hb_acc_2_dssp_index

        return

    def __repr__(self):
        def f_aa(res):
            if res == None:
                return '.'
            elif type(res) == Bio.PDB.Residue.Residue:
                return res.single_aa
            else:
                return "-"

        repr_mat = numpy.vectorize(f_aa)(numpy.array(self))

        buf = []

        for row in repr_mat:
            buf.append("".join(row))

        return os.linesep.join(buf)

    def __len__(self):
        return len(self.sheet_matrix)

    def __getitem__(self, key):
        return self.sheet_matrix[key]

    def iter_residues(self):
        for res in self.sheet_matrix.iter_residues():
            yield res

    def iter_strands(self):
        return self.sheet_matrix.iter_strands()

    def iter_strand_pairs(self):
        """Iterator over the native strand pair numbers (as tuples).

        """

        return self.sheet_matrix.iter_strand_pairs()

    def iter_bps(self, duplicates=False, offset=0):
        """

        Returns:
            iterator, where each value is of the form:
                <res_i AA>
                <res_j AA>
                <orientation i.e. "p" | "ap">
                <res_i X1 angle>
                <res_j X1 angle>
                <offset e.g. i, j + offset> 
                <res_i matrix coord>
                <res_j matrix coord>

        """

        NON_RESIDUES = set([None, "-"])

        def get_orient(res_a, res_b):
            if res_a.direction == res_b.direction:
                orient = "p"
            else:
                orient = "ap"

            return orient

        def get_hb_type(res_a, res_b):
            """

            Returns:
                a 3-tuple of (<"ap"/"p">, <"HB"/"nHB">, <"HB"/"nHB">) denoting the
                orientation and whether or not res_a and res_b are hydrogen bonded or
                not.

            """

            def get_ap_hb_type(res_a, res_b):
                """
                """

                if ((res_a.hb_don_1_dssp_index == res_a.hb_acc_1_dssp_index == res_b.dssp_index and
                     res_a.hb_don_1_energy <= DSSP_HB_THRESH and 
                     res_a.hb_acc_1_energy <= DSSP_HB_THRESH) or 
                    (res_a.hb_don_1_dssp_index == res_a.hb_acc_2_dssp_index == res_b.dssp_index and
                     res_a.hb_don_1_energy <= DSSP_HB_THRESH and 
                     res_a.hb_acc_2_energy <= DSSP_HB_THRESH) or 
                    (res_a.hb_don_2_dssp_index == res_a.hb_acc_1_dssp_index == res_b.dssp_index and
                     res_a.hb_don_2_energy <= DSSP_HB_THRESH and 
                     res_a.hb_acc_1_energy <= DSSP_HB_THRESH) or 
                    (res_a.hb_don_2_dssp_index == res_a.hb_acc_2_dssp_index == res_b.dssp_index and
                     res_a.hb_don_2_energy <= DSSP_HB_THRESH and 
                     res_a.hb_acc_2_energy <= DSSP_HB_THRESH)):
                    return ("HB", "HB")

                return ("nHB", "nHB")

            def get_p_hb_type(res_a, res_b):
                """
                """

                if ((res_a.hb_don_1_dssp_index == res_b.dssp_index - 1 and
                     res_a.hb_acc_1_dssp_index == res_b.dssp_index + 1 and
                     res_a.hb_don_1_energy <= DSSP_HB_THRESH and 
                     res_a.hb_acc_1_energy <= DSSP_HB_THRESH) or 
                    (res_a.hb_don_1_dssp_index == res_b.dssp_index - 1 and
                     res_a.hb_acc_2_dssp_index == res_b.dssp_index + 1 and
                     res_a.hb_don_1_energy <= DSSP_HB_THRESH and 
                     res_a.hb_acc_2_energy <= DSSP_HB_THRESH) or 
                    (res_a.hb_don_2_dssp_index == res_b.dssp_index - 1 and
                     res_a.hb_acc_1_dssp_index == res_b.dssp_index + 1 and
                     res_a.hb_don_2_energy <= DSSP_HB_THRESH and
                     res_a.hb_acc_1_energy <= DSSP_HB_THRESH) or
                    (res_a.hb_don_2_dssp_index == res_b.dssp_index - 1 and
                     res_a.hb_acc_2_dssp_index == res_b.dssp_index + 1 and
                     res_a.hb_don_2_energy <= DSSP_HB_THRESH and
                     res_a.hb_acc_2_energy <= DSSP_HB_THRESH)):
                    return ("HB", "nHB")

                return ("nHB", "HB")

            if res_a.direction == res_b.direction:
                pair_type = get_p_hb_type(res_a, res_b)
            else:
                pair_type = get_ap_hb_type(res_a, res_b)

            return pair_type

        def get_ca_dist(res_a, res_b):
            try:
                dist = numpy.linalg.norm(res_a["CA"].coord - res_b["CA"].coord)
            except:
                dist = "NA"

            return dist

        def get_strand_location(res):
            """

            Returns:
                True if the residue is located on an edge beta-strand i.e. it has a
                secondary structure assignment of "E" or "B" and has only on bridge
                partner.

            """

            if res.ss == 'E' or res.ss == 'B':
                if (res.bp1 and not res.bp2) or \
                   (res.bp2 and not res.bp1):

                    if res.bp1:
                        bp = res.bp1
                    else:
                        bp = res.bp2

                    return "e"
                elif res.bp1 and res.bp2:
                    return "i"
                else:
                    return "NA"

            return "NA"

        nrows = len(self)
        ncols = len(self[0])

        for i in xrange(nrows):
            if (not duplicates) and (i == nrows - 1):
                break

            for j in xrange(ncols):
                if self[i][j] == None or self[i][j] == '-':
                    continue

                for off in xrange(-offset, offset + 1):
                    if not (0 <= (j + off) < ncols):
                        continue

                    res_a = self[i][j]

                    if i > 0 and duplicates:
                        if self[i-1][j+off] not in NON_RESIDUES:
                            res_b = self[i-1][j+off]
                            orient = get_orient(res_a, res_b)
                            bridge_type = get_hb_type(res_a, res_b)
                            ca_dist = get_ca_dist(res_a, res_b)

                            yield (res_a.single_aa,
                                   res_b.single_aa,
                                   res_a.X1_angle,
                                   res_b.X1_angle,
                                   res_a.pdb_resid,
                                   res_b.pdb_resid,
                                   self.pdb_id,
                                   self.chain_id,
                                   self.sheet_id,
                                   orient,
                                   bridge_type,
                                   ca_dist,
                                   get_strand_location(res_a),
                                   get_strand_location(res_b),
                                   off,
                                   (i, j),
                                   (i - 1, j + off))

                    if i < nrows - 1 and (self[i+1][j+off] not in NON_RESIDUES):
                        res_b = self[i+1][j+off]
                        orient = get_orient(res_a, res_b)
                        bridge_type = get_hb_type(res_a, res_b)
                        ca_dist = get_ca_dist(res_a, res_b)

                        yield (res_a.single_aa,
                               res_b.single_aa,
                               res_a.X1_angle,
                               res_b.X1_angle,
                               res_a.pdb_resid,
                               res_b.pdb_resid,
                               self.pdb_id,
                               self.chain_id,
                               self.sheet_id,
                               orient,
                               bridge_type,
                               ca_dist,
                               get_strand_location(res_a),
                               get_strand_location(res_b),
                               off,
                               (i, j),
                               (i + 1, j + off))

        return

    def iter_pair_conf_tab_rows(self, offset=0):
        """

        """

        for vals in self.iter_bps(offset=offset):
            aa1, aa2, \
            X1_1, X1_2, \
            resid1, resid2, \
            pdb_id, chain_id, \
            sheet_id, \
            orient, (hb1, hb2), \
            ca_dist, loc1, loc2, \
            offset, _, _ = vals

            if X1_1 == None:
                X1_1 = "NA"
            else:
                X1_1 = "%.3f" % X1_1

            if X1_2 == None:
                X1_2 = "NA"
            else:
                X1_2 = "%.3f" % X1_2

            if ca_dist != "NA":
                ca_dist = "%.3f" % ca_dist

            row_str = "\t".join([
                aa1, aa2,
                X1_1, X1_2,
                resid1, resid2,
                pdb_id, chain_id,
                sheet_id,
                orient, hb1, hb2,
                ca_dist, loc1, loc2,
                str(offset)])

            yield row_str

        return

    def iter_bps_tab_rows(self, duplicates=False, offset=0):
        """
        """

        for bp in self.iter_bps(duplicates, offset):
            (aa1, aa2, X1, X2, orient, (aa1_hb, aa2_hb),
             ca_dist, aa1_loc, aa2_loc, off, _, _) = bp

            X1_str = "NA" if X1 == None else "%.2f" % X1
            X2_str = "NA" if X2 == None else "%.2f" % X2
            ca_dist_str = "NA" if ca_dist == None else "%.2f" % ca_dist

            yield "\t".join((aa1, aa2, X1_str, X2_str, orient, \
                    aa1_hb, aa2_hb, ca_dist_str, aa1_loc, aa2_loc, \
                    str(off)))

        return

    def get_bps_tab_str(self, duplicates=False, offset=0):
        """
        """

        return os.linesep.join(self.iter_bps_tab_rows(duplicates, offset))

    @property
    def json_dict(self):
        json_dict = {}

        json_dict["sheet_id"] = self.sheet_id
        json_dict["mat"] = self.aa_mat.tolist()
        json_dict["pdb_resid_mat"] = self.pdb_resid_mat.tolist()

        json_dict["topology"] = self.topology
        json_dict["directions"] = \
                dict((s.strand_num, s.direction) for s in self.iter_strands())
        json_dict["n_strands"] = self.n_strands
        json_dict["max_strand_length"] = self.max_strand_length

        json_dict["is_native"] = True
        json_dict["is_bifurcated"] = self.is_bifurcated
        json_dict["is_barrel"] = self.is_barrel
        json_dict["is_bulged"] = self.is_bulged

        strands_dict = {}

        for strand_num, strand in self.sheet_matrix.sheet.strand_ord_dict.iteritems():
            strand_str = " ".join(res.global_id for res in strand)
            strands_dict[strand.strand_num] = strand_str

        json_dict["strands"] = strands_dict

        residues_dict = {}

        for res in self.iter_residues():
            res_dict = { "name" : res.single_aa,
                         "direction" : res.direction }
            residues_dict[res.pdb_resid] = res_dict

        json_dict["residues"] = residues_dict

        return json_dict

    @property
    def json_string(self):
        return json.dumps(self.json_dict)

    def get_strand_pair_json_dict(self, strand_pair):
        """Get a JSON string representation of a beta-matrix
        subset containing only the strands specified in
        <strand_pair>.

        """

        s1, s2 = strand_pair

        json_dict = self.json_dict

        json_dict["n_strands"] = 2
        json_dict["max_strand_length"] = max(len(self.strands[s1]),
                                             len(self.strands[s2]))
        json_dict["is_barrel"] = False
        json_dict["is_native"] = True

        directions = json_dict["directions"]

        json_dict["directions"] = \
                dict( (s, d) for (s, d) in directions.iteritems()
                      if s in strand_pair )

        strands_dict = json_dict["strands"]

        json_dict["strands"] = \
                dict( (s, ss) for (s, ss) in strands_dict.iteritems()
                      if s in strand_pair )

        residues = {}

        for strand in self.strands.itervalues():
            if strand.strand_num not in strand_pair:
                continue

            for res in strand:
                residues[res.pdb_resid] = { "name" : res.single_aa,
                                            "direction" : res.direction }

        new_topology = []

        for group in self.topology:
            group_lst = []

            for strand_num in group:
                if strand_num not in strand_pair:
                    continue

                group_lst.append(strand_num)

            if group_lst:
                new_topology.append(group_lst)

        aa_mat = []
        pdb_resid_mat = []

        seen_pdb_resids = set()

        # Generate the plain-text matrix.
        for i, row in enumerate(self):
            aa_row_list = []
            pdb_resid_row_list = []
            row_has_aas = False
            row_is_barrel_edge = False

            for j, val in enumerate(row):
                if val == None:
                    aa_row_list.append( "." )
                    pdb_resid_row_list.append( "." )
                elif val == '-':
                    aa_row_list.append( "-" )
                    pdb_resid_row_list.append( "-" )
                else:
                    aa = self[i][j].single_aa
                    pdb_resid = self[i][j].pdb_resid

                    if pdb_resid in seen_pdb_resids:
                        continue

                    seen_pdb_resids.add(pdb_resid)

                    if self[i][j].strand_num not in strand_pair:
                        aa_row_list.append(".")
                        pdb_resid_row_list.append(".")
                        continue

                    row_has_aas = True
                    aa_row_list.append(aa)
                    pdb_resid_row_list.append(pdb_resid)

            if not row_has_aas:
                continue

            aa_mat.append(aa_row_list)
            pdb_resid_mat.append(pdb_resid_row_list)

        _aa_mat = numpy.array(aa_mat)
        _pdb_resid_mat = numpy.array(pdb_resid_mat)

        # Remove empty columns in the matrix.
        cols_to_delete = []

        for col in xrange(_aa_mat.shape[1]):
            if numpy.all(_aa_mat[:,col] == ['.', '.']):
                cols_to_delete.append(col)

        aa_mat = numpy.delete(_aa_mat, cols_to_delete, 1)
        pdb_resid_mat = numpy.delete(_pdb_resid_mat, cols_to_delete, 1)

        json_dict["sheet_id"] += "-%02d-%02d" % strand_pair
        json_dict["mat"] = aa_mat.tolist()
        json_dict["pdb_resid_mat"] = pdb_resid_mat.tolist()
        json_dict["is_bulged"] = True if "-" in aa_mat else False
        json_dict["is_bifurcated"] = False
        json_dict["topology"] = new_topology
        json_dict["residues"] = residues

        return json_dict

    def iter_strand_pair_json_strings(self):
        """
        """

        for pair in self.iter_strand_pairs():
            json_dict = self.get_strand_pair_json_dict(pair)
            yield json.dumps(json_dict)

        return
