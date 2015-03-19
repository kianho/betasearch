#!/usr/bin/env python
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Description:
    ...

"""

import os
import sys
import numpy
import networkx as nx
import cPickle
import pprint
import copy
import re
import time

from collections import defaultdict, deque
from itertools import product

NO_OVERLAP = 0
BRIDGE = 1
PEPTIDE = 2
OPP_BRIDGE = 4
OPP_PEPTIDE = 5

L_TRIMER = 1
V_TRIMER = 3
H_TRIMER = 5
SYM_V_TRIMER = 15
SYM_H_TRIMER = 31

WC_RE = re.compile("[\*]")

AAS = "ARNDCQEGHILKMFPSTWYV"


class Elbow:
    """This class stores the row and column indices of a trimer "elbow", which is
    the central residue in the trimer.
    
    For example, 'B' is the elbow of the L-trimer below:

       0 1 2            
      0. A .
      1. B C
         ^
         :.. Elbow @ row 1, column 1

    """

    def __init__(self, row, col):
        """Constructor.

        Parameters
        ----------
        row : int
            Row co-ordinate of the elbow.
        col : int
            Column co-ordinate of the elbow.

        """
        self.row = row
        self.col = col

    def same(self, other):
        """Check if this elbow is in the same position as the other.  This
        assumes that the parent trimers of both elbows are from the same
        beta-matrix.

        """
        return self.row == other.row and self.col == other.col


class Span:
    """This is a glorified tuple denoting the row or column extents of a given
    trimer. It is typically used to determine if two trimers overlap in a
    beta-matrix.

    """

    def __init__(self, begin, end):
        """Constructor.

        Parameters
        ----------
        begin : int
            Typically, this is the row/column index of the elbow.
        end : int
            The outer row/column index of the trimer.

        """

        self.begin = begin
        self.end = end

    def get_min(self):
        if self.begin < self.end:
            return self.begin
        return self.end

    def get_max(self):
        if self.begin > self.end:
            return self.begin
        return self.end

    def get_reversed(self):
        return Span(self.end, self.begin)

    def get_tuple(self):
        return (self.begin, self.end)

    def same_direction(self, other):
        """Check if the "direction" of two  spans are the same. For L-trimers,
        directions beginning from the elbow outwards, horizontally
        (for column spans) or vertically for (row spans).

        For example:

                    trimer A        trimer B

                                    +ve column span direction 
                                    <-- .
             +ve  ^ A               C B | -ve row span direction
             row  | B C               A v
             span . -->
             dir.   -ve column span direction

        Therefore, the corresponding row spans and column spans are in opposite
        directions between both trimers.

        """
        return ((self.begin - self.end) < 0) == ((other.begin - other.end) < 0)

    def overlaps(self, other):
        if self.begin == other.begin and self.end == other.end:
            return True
        if self.begin == other.end and self.end == other.begin:
            return True

        this_begin = self.get_min()
        this_end = self.get_max()

        other_begin = other.get_min()
        other_end = other.get_max()

        return (((this_begin < other_begin) and (other_begin < this_end)) or
                ((other_begin < this_begin) and (this_begin < other_end)))


class Trimer:
    """
    """

    def __init__(self, type_, seq, orient, elbow, span1, span2, num):
        """Constructor.

        """

        self.type_ = type_
        self.seq = seq
        self.orient = orient
        self.elbow = elbow
        self.span1 = span1
        self.span2 = span2
        self.num = num

        self.id_str = "%d_%s" % (self.type_, seq)

        if type_ == L_TRIMER:
            self.equiv_orients = 1 << orient
        elif type_ == V_TRIMER:
            if orient & 2:
                self.equiv_orients = 12
            else:
                self.equiv_orients = 3
        elif type_ == H_TRIMER:
            if orient & 1:
                self.equiv_orients = 10
            else:
                self.equiv_orients = 5
        else:
            self.equiv_orients = 15

        return

    def __repr__(self):
        return "%s:%d" % (self.id_str, self.num)

    def __hash__(self):
        return self.num

    def __eq__(self, other):
        return  self.__hash__() == hash(other) 

    def get_L_overlap_type(self, other):
        if other.type_ == L_TRIMER:
            if self.elbow.same(other.elbow):
                rel_orient = self.orient ^ other.orient

                if rel_orient == 2:
                    return PEPTIDE
                elif rel_orient == 1:
                    return BRIDGE
                else:
                    return NO_OVERLAP
            elif self.elbow.row == other.elbow.row:
                if self.span2.overlaps(other.span2):
                    return OPP_PEPTIDE
            elif self.elbow.col == other.elbow.col:
                if self.span1.overlaps(other.span1):
                    return OPP_BRIDGE

        elif other.type_ == V_TRIMER or other.type_ == SYM_V_TRIMER:
            if self.elbow.col == other.elbow.col:
                if self.span1.overlaps(other.span1):
                    other_span = other.span1
                elif self.span1.overlaps(other.span2):
                    other_span = other.span2
                else:
                    return NO_OVERLAP

                if self.span1.same_direction(other_span):
                    return BRIDGE
                else:
                    return OPP_BRIDGE
        else:
            if self.elbow.row == other.elbow.row:
                if self.span2.overlaps(other.span1):
                    other_span = other.span1
                elif self.span2.overlaps(other.span2):
                    other_span = other.span2
                else:
                    return NO_OVERLAP

                if self.span2.same_direction(other_span):
                    return PEPTIDE
                else:
                    return OPP_PEPTIDE

        return NO_OVERLAP

    def get_V_overlap_type(self, other):
        if self.elbow.col != other.elbow.col:
            return NO_OVERLAP

        if other.type_ == V_TRIMER or other.type_ == SYM_V_TRIMER:
            if self.span1.overlaps(other.span1):
                this_span = self.span1
                other_span = other.span1
            elif self.span1.overlaps(other.span2):
                this_span = self.span1
                other_span = other.span2
            elif self.span2.overlaps(other.span1):
                this_span = self.span2
                other_span = other.span1
            elif self.span2.overlaps(other.span2):
                this_span = self.span2
                other_span = other.span2
            else:
                return NO_OVERLAP

            if this_span.same_direction(other_span):
                return BRIDGE
            else:
                return OPP_BRIDGE
        elif other.type_ == L_TRIMER:
            if self.span1.overlaps(other.span1):
                this_span = self.span1
                other_span = other.span1
            elif self.span2.overlaps(other.span1):
                this_span = self.span2
                other_span = other.span1
            else:
                return NO_OVERLAP

            if this_span.same_direction(other_span):
                return BRIDGE
            else:
                return OPP_BRIDGE

        return NO_OVERLAP

    def get_H_overlap_type(self, other):
        if self.elbow.row != other.elbow.row:
            return NO_OVERLAP

        if other.type_ == H_TRIMER or other.type_ == SYM_H_TRIMER:
            if self.span1.overlaps(other.span1):
                this_span = self.span1 
                other_span = other.span1 
            elif self.span1.overlaps(other.span2):
                this_span = self.span1 
                other_span = other.span2 
            elif self.span2.overlaps(other.span1):
                this_span = self.span2
                other_span = other.span1
            elif self.span2.overlaps(other.span2):
                this_span = self.span2
                other_span = other.span2
            else:
                return NO_OVERLAP

            if this_span.same_direction(other_span):
                return PEPTIDE
            else:
                return OPP_PEPTIDE

        elif other.type_ == L_TRIMER:
            if self.span1.overlaps(other.span2):
                this_span = self.span1
                other_span = other.span2
            elif self.span2.overlaps(other.span2):
                this_span = self.span2
                other_span = other.span2
            else:
                return NO_OVERLAP
           
            if this_span.same_direction(other_span):
                return PEPTIDE
            else:
                return OPP_PEPTIDE

        return NO_OVERLAP

    def get_overlap_type(self, other):
        if self.type_ == 1:
            return self.get_L_overlap_type(other)
        elif self.type_ == 3 or self.type_ == 15:
            return self.get_V_overlap_type(other)
        elif self.type_ == 5 or self.type_ == 31:
            return self.get_H_overlap_type(other)

        return NO_OVERLAP

    def overlapping_span_num(self, other):
        overlap = self.get_overlap_type(other)

        assert(overlap != NO_OVERLAP)

        if self.type_ == L_TRIMER:
            if overlap == BRIDGE or overlap == OPP_BRIDGE:
                return 1
            else:
                return 2
        elif self.type_ == V_TRIMER or self.type_ == SYM_V_TRIMER:
            if other.type_ == L_TRIMER:
                if self.span1.overlaps(other.span1):
                    return 1
                elif self.span2.overlaps(other.span1):
                    return 2
                else:
                    return 0
            else:
                assert(other.type_ == V_TRIMER or
                       other.type_ == SYM_V_TRIMER)

                # other trimer _must_ be a V-trimer.
                if self.span1.overlaps(other.span1) or \
                   self.span1.overlaps(other.span2):
                    return 1
                elif self.span2.overlaps(other.span1) or \
                     self.span2.overlaps(other.span2):
                    return 2
                else:
                    return 0
        elif self.type_ == H_TRIMER or self.type_ == SYM_H_TRIMER:
            if other.type_ == L_TRIMER:
                if self.span1.overlaps(other.span2):
                    return 1 
                elif self.span2.overlaps(other.span2):
                    return 2 
                else:
                    return 0 
            else:
                assert(other.type_ == H_TRIMER or
                       other.type_ == SYM_H_TRIMER)

                # other trimer _must_ be an H-trimer.
                if self.span1.overlaps(other.span1) or \
                   self.span1.overlaps(other.span2):
                    return 1
                elif self.span2.overlaps(other.span1) or \
                     self.span2.overlaps(other.span2):
                    return 2
                else:
                    return 0 

        return 0


class Query:
    def __init__(self, line):
        """Constructor.

        """

        if len(line.split(":")) == 1:
          self._line = line
        else:
          self._line = line.split(":")[-1]

        self._line = line.strip()
        self._first_l_trimer = None
        self._first_asym_trimer = None

        self._init_matrix()
        self._build_trimer_structs()
        self._build_overlaps()

        return

    def _init_matrix(self):
        self._mat = numpy.array(map(lambda row : list(row),
                                self._line.split(":")[-1].strip().split(",")))
        return

    def _build_trimer_structs(self):
        """
        """

        trimers = []

        bridges = defaultdict(lambda : defaultdict(list))
        peptides = defaultdict(lambda : defaultdict(list))

        for trimer in self.trimers_gen():
            if self._first_l_trimer == None:
                if trimer.type_ == L_TRIMER:
                    self._first_l_trimer = trimer

            if self._first_asym_trimer == None:
                if trimer.type_ % 16 != 15:
                    self._first_asym_trimer = trimer

            trimers.append(trimer)

            begin1, end1 = trimer.span1.get_min(), trimer.span1.get_max()
            begin2, end2 = trimer.span2.get_min(), trimer.span2.get_max()

            if trimer.type_ == 1:
                bridges[trimer.elbow.col][(begin1, end1)].append(trimer)
                peptides[trimer.elbow.row][(begin2, end2)].append(trimer)
            elif trimer.type_ == 3 or trimer.type_ == 15:
                bridges[trimer.elbow.col][(begin1, end1)].append(trimer)
                bridges[trimer.elbow.col][(begin2, end2)].append(trimer)
            else:
                peptides[trimer.elbow.row][(begin1, end1)].append(trimer)
                peptides[trimer.elbow.row][(begin2, end2)].append(trimer)

        self._bridges = bridges
        self._peptides = peptides
        self._trimers = trimers

        return

    def _build_overlaps(self):
        """Construct the trimer overlap graph.

        """

        G = nx.DiGraph()

        for val in self._bridges.itervalues():
            for trimers in val.itervalues():
                for i in xrange(len(trimers)):
                    for j in xrange(i + 1, len(trimers)):
                        src, des = trimers[i], trimers[j]

                        G.add_edge(src, des,
                                   span_num=src.overlapping_span_num(des),
                                   rel_orient=src.orient ^ des.orient,
                                   overlap_type=src.get_overlap_type(des))

                        G.add_edge(des, src,
                                   span_num=des.overlapping_span_num(src),
                                   rel_orient=des.orient ^ src.orient,
                                   overlap_type=des.get_overlap_type(src))


        for val in self._peptides.itervalues():
            for trimers in val.itervalues():
                for i in xrange(len(trimers)):
                    for j in xrange(i + 1, len(trimers)):
                        src, des = trimers[i], trimers[j]

                        G.add_edge(src, des,
                                   span_num=src.overlapping_span_num(des),
                                   rel_orient=src.orient ^ des.orient,
                                   overlap_type=src.get_overlap_type(des))

                        G.add_edge(des, src,
                                   span_num=des.overlapping_span_num(src),
                                   rel_orient=des.orient ^ src.orient,
                                   overlap_type=des.get_overlap_type(src))

        self._overlaps = G

        return

    def get_whoosh_query_str(self, disjunctive=False):
        if disjunctive:
            operator = " OR "
        else:
            operator = " AND "

        return unicode(operator.join(t.id_str for t in self._trimers))

    def get_pretty_str(self):
        return line_to_pretty_str(self._line)

    def trimers_gen(self):
        """Generator over all the trimers in a beta-matrix.

        """

        mat = self._mat
        num = 0
        num_rows, num_cols = mat.shape

        for row in xrange(num_rows):
            for col in xrange(num_cols):
                if row < num_rows - 1 and col < num_cols - 1:
                    for orient in xrange(4):
                        trimer = make_L_trimer(row, col, orient, mat, num)

                        if trimer != None:
                            num += 1
                            yield trimer

                if row < num_rows - 2:
                    trimer = make_V_trimer(row, col, mat, num)

                    if trimer != None:
                        num += 1
                        yield trimer

                if col < num_cols - 2:
                    trimer = make_H_trimer(row, col, mat, num)

                    if trimer != None:
                        num += 1
                        yield trimer

        return

    def verify(self, results, trimers_dir):
        if len(results) == 0:
            return

        # If there is only one trimer in the query, then there is no need to
        # filter/verify the results.
        if len(self._trimers) == 1:
            for hit in results:
                yield (hit["sheet_id"], hit["molecule_name"],
                        hit["organism_common_name"],
                        hit["organism_scientific_name"])
                #yield str(hit["sheet_id"])
            return

        if self._first_l_trimer != None:
            root = self._first_l_trimer
        elif self._first_asym_trimer != None:
            root = self._first_asym_trimer
        else:
            root = self._trimers[0]

        for hit in results:
            sheet_id = str(hit["sheet_id"]) 
            molecule_name = str(hit["molecule_name"])
            organism_common_name = str(hit["organism_common_name"])
            organism_scientific_name = str(hit["organism_scientific_name"])
            found = False

            f = open(os.path.join(trimers_dir, "%s.record" % sheet_id))
            disk_record = cPickle.load(f)
            f.close()

            for root_src in disk_record["trimers"][root.id_str]:
                if root_src.type_ % 16 == 15:
                    max_span_combos = 2
                else:
                    max_span_combos = 1

                # For each span combination.
                for combo in xrange(max_span_combos):
                    trimer_spans = {}

                    matched_trimers = deque(maxlen=len(self._trimers))
                    matched_trimers.append((root, root_src))

                    if combo == 0:
                        trimer_spans[root_src] = { 1 : root_src.span1,
                                                   2 : root_src.span2 }
                    else:
                        trimer_spans[root_src] = { 2 : root_src.span1,
                                                   1 : root_src.span2 }

                    visited = defaultdict(bool)
                    t_visited = defaultdict(bool)

                    q_des_to_t_src = {}
                    q_des_to_q_src = {}
                    queue = deque(maxlen=len(self._trimers))

                    # Push the root's neighbours onto the queue.
                    for neigh in self._overlaps.neighbors_iter(root):
                        q_des_to_t_src[neigh] = root_src
                        q_des_to_q_src[neigh] = root
                        queue.append(neigh)

                    q_src = root
                    visited[q_src] = True
                    found = True

                    while len(queue) > 0:
                        q_des = queue.popleft()
                        t_src = q_des_to_t_src[q_des]

                        t_visited[t_src] = True

                        assert(t_src != None)

                        if visited[q_des]:
                            continue

                        if not self._overlaps.has_edge(q_src, q_des):
                            q_src = q_des_to_q_src[q_des]

                        visited[q_des] = True

                        src_span_num = self._overlaps[q_src][q_des]["span_num"]
                        des_span_num = self._overlaps[q_des][q_src]["span_num"]

                        rel_orient = self._overlaps[q_src][q_des]["rel_orient"]
                        des_orient = t_src.orient ^ rel_orient

                        equiv_orients = get_equiv_orients(des_orient,
                                                          q_des.type_)

                        src_span = trimer_spans[t_src][src_span_num]
                        t_des_id_str = q_des.id_str
                        overlap_type = \
                            self._overlaps[q_src][q_des]["overlap_type"]

                        if overlap_type == OPP_PEPTIDE:
                            des_span = src_span.get_reversed()
                            index = disk_record["col-index"]
                            rowcol = t_src.elbow.row
                        elif overlap_type == PEPTIDE:
                            des_span = src_span
                            index = disk_record["col-index"]
                            rowcol = t_src.elbow.row
                        elif overlap_type == OPP_BRIDGE:
                            des_span = src_span.get_reversed()
                            index = disk_record["row-index"]
                            rowcol = t_src.elbow.col
                        elif overlap_type == BRIDGE:
                            des_span = src_span
                            index = disk_record["row-index"]
                            rowcol = t_src.elbow.col
                        else:
                            pass

                        t_des = index.get((rowcol, des_span.get_tuple(), t_des_id_str,
                                           equiv_orients), None)

                        if t_des == None:
                            found = False
                            break

                        if t_visited[t_des]:
                            found = False
                            break
                        else:
                            t_visited[t_des] = True

                        t_des.orient = des_orient
                        trimer_spans[t_des] = { des_span_num : des_span }

                        if (des_span_num == 1):
                            trimer_spans[t_des][2] = copy.deepcopy(t_des.span2)
                        else:
                            trimer_spans[t_des][1] = copy.deepcopy(t_des.span1)

                        for neigh in self._overlaps.neighbors_iter(q_des):
                            if visited[neigh]:
                                continue

                            q_des_to_t_src[neigh] = t_des
                            q_des_to_q_src[neigh] = q_des
                            queue.append(neigh)

                        matched_trimers.append((q_des, t_des))

                        q_src = q_des

                    if found:
                        break

                if found:
                    break

            if found:
                yield sheet_id, molecule_name, organism_common_name, organism_scientific_name
#               t1 = time.clock()
#               yield sheet_id, t1 - t0
#               t0 = time.clock()
                
        return


def get_equiv_orients(orient, trimer_type):
    """
    """

    if trimer_type == L_TRIMER:
        return 1 << orient
    elif trimer_type == V_TRIMER:
        if orient & 2:
            return 12
        else:
            return 3
    elif trimer_type == H_TRIMER:
        if orient & 1:
            return 10
        else:
            return 5

    return 15


def get_neighbour_col(row, col, direction, mat):
    new_col = col

    try:
        while mat[row][new_col] == '-':
            new_col += direction
    except IndexError:
        return None

    return new_col


def make_L_trimer(row, col, orient, mat, num):
    """
    """

    if orient == 0:
        elbow = Elbow(row, col)
        span1 = Span(elbow.row, elbow.row + 1)
        span2 = Span(elbow.col, elbow.col + 1)
    elif orient == 1:
        elbow = Elbow(row, col + 1)
        span1 = Span(elbow.row, elbow.row + 1)
        span2 = Span(elbow.col, elbow.col - 1)
    elif orient == 2:
        elbow = Elbow(row + 1, col)
        span1 = Span(elbow.row, elbow.row - 1)
        span2 = Span(elbow.col, elbow.col + 1)
    elif orient == 3:
        elbow = Elbow(row + 1, col + 1)
        span1 = Span(elbow.row, elbow.row - 1)
        span2 = Span(elbow.col, elbow.col - 1)

    seq = [mat[span1.end][elbow.col],
           mat[elbow.row][elbow.col],
           mat[elbow.row][span2.end]]

    if seq[0] == '.' or seq[0] == '-' or \
       seq[1] == '.' or seq[1] == '-' or \
       seq[2] == '.':
        return

    if seq[2] == '-':
        if orient == 1 or orient == 3:
            direction = -1
        else:
            direction = 1

        span2.end = get_neighbour_col(elbow.row, span2.end, direction, mat) 
        seq[2] = mat[elbow.row][span2.end]
    
    return Trimer(1, "".join(seq), orient, elbow, span1, span2, num)


def make_V_trimer(row, col, mat, num):
    """
    """

    seq = [None, None, None]

    for i in xrange(3):
        seq[i] = mat[row + i][col]

        if seq[i] == '.' or seq[i] == '-':
            return

    elbow = Elbow(row + 1, col)
    span1 = Span(elbow.row, elbow.row - 1)
    span2 = Span(elbow.row, elbow.row + 1)

    type_ = 3

    if seq[0] < seq[2]:
        orient = 0
    elif seq[0] > seq[2]:
        orient = 2
        seq.reverse()
        begin, end = span1.begin, span1.end
        span1.begin, span1.end = span2.begin, span2.end
        span2.begin, span2.end = begin, end
    else:
        orient = 0
        type_ = 15

    return Trimer(type_, "".join(seq), orient, elbow, span1, span2, num)


def make_H_trimer(row, col, mat, num):
    """
    """

    elbow = Elbow(row, col + 1)
    span1 = Span(elbow.col, elbow.col - 1)
    span2 = Span(elbow.col, elbow.col + 1)

    seq = [mat[elbow.row][elbow.col - 1],
           mat[elbow.row][elbow.col],
           mat[elbow.row][elbow.col + 1]]

    if seq[0] == '.' or seq[1] == '.' or \
       seq[1] == '-' or seq[2] == '.':
        return

    if seq[0] == '-':
        span1.end = get_neighbour_col(row, span1.end, -1, mat)
    
    if seq[2] == '-':
        span2.end = get_neighbour_col(row, span2.end, 1, mat)

    seq[0] = mat[row][span1.end]
    seq[2] = mat[row][span2.end]

    type_ = 5

    if seq[0] < seq[2]:
        orient = 0
    elif seq[0] > seq[2]:
        seq.reverse()
        begin, end = span1.begin, span1.end
        span1.begin, span1.end = span2.begin, span2.end
        span2.begin, span2.end = begin, end
        orient = 1
    else:
        orient = 0
        type_ = 31

    return Trimer(type_, "".join(seq), orient, elbow, span1, span2, num)


def update_disk_record(record, trimer):
    """
    """

    if trimer.type_ == L_TRIMER:
        r_span = (trimer.span1.begin, trimer.span1.end)
        c_span = (trimer.span2.begin, trimer.span2.end)

        key1 = (trimer.elbow.col, r_span, trimer.id_str, trimer.equiv_orients)
        key2 = (trimer.elbow.row, c_span, trimer.id_str, trimer.equiv_orients)

        record["row-index"][key1] = trimer
        record["col-index"][key2] = trimer

    elif trimer.type_ == V_TRIMER or trimer.type_ == SYM_V_TRIMER:
        r_span1 = (trimer.span1.begin, trimer.span1.end)
        r_span2 = (trimer.span2.begin, trimer.span2.end)

        key1 = (trimer.elbow.col, r_span1, trimer.id_str, trimer.equiv_orients)
        key2 = (trimer.elbow.col, r_span2, trimer.id_str, trimer.equiv_orients)

        record["row-index"][key1] = trimer
        record["row-index"][key2] = trimer

    else:
        c_span1 = (trimer.span1.begin, trimer.span1.end)
        c_span2 = (trimer.span2.begin, trimer.span2.end)

        key1 = (trimer.elbow.row, c_span1, trimer.id_str, trimer.equiv_orients)
        key2 = (trimer.elbow.row, c_span2, trimer.id_str, trimer.equiv_orients)

        record["col-index"][key1] = trimer
        record["col-index"][key2] = trimer

    record["trimers"][trimer.id_str].add(trimer)
    record["trimers-list"].append(trimer)

    return


def line_to_matrix(line):
    line = line.strip()

    if len(line.split(":")) > 1:
        line = line.split(":")[-1]

    return numpy.array(map(lambda row : list(row),
                       line.strip().split(",")))


def line_to_pretty_str(line):
    """
    """

    line = line.strip()

    if len(line.split(":")) > 1:
        line = line.split(":")[-1]

    return os.linesep.join(row for row in line.split(','))


def line_to_num_res(line):
    return int(line.split(":")[-2])


def get_sheet_id(line):
    return ":".join(line.strip().split(":")[:2])


def is_valid_matrix(mat):
    print mat.shape


def trimers_gen(mat):
    """Generator over all the trimers in a beta-matrix.

    """

    num = 0
    num_rows, num_cols = mat.shape

    for row in xrange(num_rows):
        for col in xrange(num_cols):
            if row < num_rows - 1 and col < num_cols - 1:
                for orient in xrange(4):
                    trimer = make_L_trimer(row, col, orient, mat, num)

                    if trimer != None:
                        num += 1
                        yield trimer

            if row < num_rows - 2:
                trimer = make_V_trimer(row, col, mat, num)

                if trimer != None:
                    num += 1
                    yield trimer

            if col < num_cols - 2:
                trimer = make_H_trimer(row, col, mat, num)

                if trimer != None:
                    num += 1
                    yield trimer

    return


def make_graph(line, verbose=False):
    """

    Arguments:

    """

    g = nx.Graph()
    mat = line_to_matrix(line)
    num_rows, num_cols = mat.shape

    for row in xrange(num_rows):
        for col in xrange(num_cols):
            if mat[row, col] == "." or mat[row, col] == "-":
                continue
            g.add_node((row, col), aa=mat[row, col])

    for row in xrange(num_rows):
        for col in xrange(num_cols):
            if re.match(r"[.-]", mat[row,col]):
                continue

            # add right peptide edge
            if col < num_cols - 1:
                next_col = get_neighbour_col(row, col + 1, 1, mat)

                if mat[row, next_col] != '.':
                    src = (row, col)
                    des = (row, next_col)
                    g.add_edge(src, des)

            # add bottom bridge edge
            if row < num_rows - 1 and not re.match(r"[.-]", mat[row + 1][col]):
                src = (row, col)
                des = (row + 1, col)
                g.add_edge(src, des)

    if verbose:
        pprint.pprint([g.node[n]["aa"] for n in g.nodes_iter()])
        pprint.pprint([(g.node[src]["aa"], g.node[des]["aa"])
                       for (src, des) in g.edges_iter()])

    return g, mat


def line_is_connected(line):
    """

    Returns:
        True if the implied graph is connected (ie. has a path from any node to
        any other node). Otherwise, False.

    """

    g, _ = make_graph(line)
    root = g.nodes()[0]
    stack = deque([root], maxlen=g.number_of_nodes())
    visited = defaultdict(bool)

    while len(stack) > 0:
        curr = stack.popleft()
        visited[curr] = True

        for n in g.neighbors_iter(curr):
            if not visited[n]:
                stack.appendleft(n)

    # check to see if all nodes were visited.
    for v in g.nodes_iter():
        if not visited[v]:
            return False

    return True


def graph_is_connected(g):
    """

    Returns:
        True if the implied graph is connected (ie. has a path from any node to
        any other node). Otherwise, False.

    """

    root = g.nodes()[0]
    stack = deque([root], maxlen=g.number_of_nodes())
    visited = defaultdict(bool)

    while len(stack) > 0:
        curr = stack.popleft()
        visited[curr] = True

        for n in g.neighbors_iter(curr):
            if not visited[n]:
                stack.appendleft(n)

    # check to see if all nodes were visited.
    for v in g.nodes_iter():
        if not visited[v]:
            return False

    return True


def wc_expansions_iter(line):
    """This is the "smart" expansion.

    """

    wc_indices = (m.start() for m in re.finditer(WC_RE, line))

    queue = deque()
    queue.append(line)

    while len(queue) > 0:
        curr_line = queue.popleft()

        # lazy expansion.
        lazy_line, n_wcs = re.subn(WC_RE, '.', curr_line)

        if line_is_connected(lazy_line):
            # don't need to manually expand!
            yield lazy_line
        else:
            for aa in AAS:
                new_line = re.sub(WC_RE, aa, curr_line, count=1) 
                queue.append(new_line)

    return


def enforced_wc_expansions_iter(line):
    """
    """

    wc_indices = [m.start() for m in re.finditer(WC_RE, line)]

    if len(wc_indices) == 1:
        for aa in AAS:
            yield line.replace("*", aa)
    elif len(wc_indices) == 2:
        expansion_aas = (x for x in product(AAS, AAS))

        for aa1, aa2 in expansion_aas:
            yield line.replace("*", aa1, 1).replace("*", aa2, 1)
    else:
        yield line

    return
