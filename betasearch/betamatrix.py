#!/usr/bin/env python
"""

Author: Kian Ho <hui.kian.ho@gmail.com>
Description:
    This is a module containing miscellaneous operations to be performed on
    sheet matrix files (i.e. files with the *.sm extension).

"""

import os
import sys
import optparse
import copy
import itertools
import math
import glob

from collections import defaultdict, deque
from pprint import pprint
from _betasearch import *


three_to_one = defaultdict(lambda : None,
    { "ALA" : "A", "CYS" : "C", "ASP" : "D", "GLU" : "E",
      "PHE" : "F", "GLY" : "G", "HIS" : "H", "ILE" : "I",
      "LYS" : "K", "LEU" : "L", "MET" : "M", "ASN" : "N",
      "PRO" : "P", "GLN" : "Q", "ARG" : "R", "SER" : "S",
      "THR" : "T", "VAL" : "V", "TRP" : "W", "TYR" : "Y" })


class Translations:
    def __init__(self, classes_fn):
        """Constructor.

        """

        self.classes = {}
        self.class_names = {}
        self.res_to_classes = defaultdict(list)

        self.__load_classes(classes_fn)

        return

    def __load_classes(self, classes_fn):
        """
        """

        f = open(classes_fn, "r")

        for line in f:
            if line.startswith("#"):
                continue

            vals = line.split()
            class_code = vals[0]
            class_name = vals[1]

            self.classes[class_code] = \
                    map(lambda r : three_to_one[r.upper()], vals[2:])
            self.class_names[class_code] = class_name

            for res in vals[2:]:
                self.res_to_classes[three_to_one[res.upper()]].append(class_code)

        f.close()

        return


class BetaMatrix:
    """This class represents a beta-matrix that is read from a *.bm text file.
    This class will probably be merged with the BetaMatrix class in the betapy
    module at some stage.

    """

    def __init__(self, bm_fn=None, line=None):
        """Constructor.

        Arguments:
            bm_fn -- path to the beta-matrix *.bm file.

        Returns:
            None

        """

        self.bm_fn = bm_fn
        self.matrix = []
        self.header = None
        self.root_sheet_id = None
        self.sheet_id = None
        self.nrows = None
        self.ncols = None
        self.nentries = None
        self.nmutations = 0
        self.is_barrel = 'f'

        self._load_matrix(bm_fn, line)

        return

    def __repr__(self):
        buf = self.get_header() + os.linesep
        
        for row in self.matrix:
            buf += "".join(row) + os.linesep

        return buf[:-1]

    def _load_matrix(self, bm_fn=None, line=None):
        """Load the beta-matrix from the *.bm file into a two-dimensional
        list-of-lists, representing the "raw" beta-matrix that is wrapped by
        this class.

        Returns:
            None

        """

        if bm_fn != None:
            # Load the beta-matrix from a *.bm file.
            f = open(self.bm_fn, "r")
            line = f.readlines()[0].strip()
            header = parse_beta_matrix_header(line)
            self.sheet_id = header["sheet_id"]
            self.root_sheet_id = header["root_sheet_id"]
            self.matrix.extend(header["lines"])
            f.close()
        elif line != None:
            # Load the beta-matrix from a list of beta-matrix lines (assumed to
            # be read from a *.bm file.
            header = parse_beta_matrix_header(line)
            self.sheet_id = header["sheet_id"]
            self.root_sheet_id = header["root_sheet_id"]
            self.matrix.extend(header["lines"])

        self.nrows = header["n_rows"]
        self.ncols = header["n_cols"]
        self.nentries = header["n_entries"]
        self.is_barrel = header["is_barrel"]

        return

    def calc_nentries(self):
        """Calculate and set the number of entries in this BetaMatrix.

        """

        self.nentries = len(list(self.gen_entries()))

        return self.nentries

    def set(self, row, col, val):
        self.matrix[row][col] = val
        return

    def get(self, row, col):
        return self.matrix[row][col]

    def get_root_sheet_id(self):
        return self.sheet_id.split(":")[0]

    def get_sheet_id(self):
        return self.sheet_id

    def get_repr(self, style="ugly"):
        if style == "ugly":
            buf = "%s:%s" % (self.get_header(), 
                             ",".join("".join(row) for row in self.matrix))
        else:
            buf = self.__repr__()

        return buf

    def get_header(self, root_sheet_id=False):
        return "%s:%d:%d:%s:%d" % (self.get_sheet_id(), self.nrows,
                                    self.ncols, self.is_barrel, self.nentries)

    def set_sheet_id(self, sheet_id):
        self.sheet_id = sheet_id
        self.root_sheet_id = "_".join(self.sheet_id.split("_")[:3])
        return

    def clear(self):
        """Fill the entire beta-matrix with padding characters i.e. ".".

        Returns:
            None

        """

        for row in xrange(self.nrows):
            for col in xrange(self.ncols):
                self.set(row, col, ".")

        return

    def is_padding(self, row, col):
        return self.get(row, col) == "."

    def is_stretch(self, row, col):
        return self.get(row, col) == "-"

    def get_center_entry(self):
        """Get the (<row>, <col>, <amino acid>) entry located in the approximate
        middle of the beta-matrix. For example:


            matrix =

               0    1    2 
            [["A", "B", "C"], 0
             ["D", "E", "F"], 1
             ["G", "H", "I"]] 2

        the center entry is located at matrix[1][1], however, for beta-matrices
        that do not have an odd-number of rows and cols, the center entry is
        approximated as follows:

            row = floor(nrows / 2)
            col = floor(ncols / 2)

        therefore if

            matrix = 

               0    1
            [["A", "B"], 0
             ["C", "D"]] 1

        then the center will be located at matrix[1][1].

        Returns:
            a tuple of the form (<row>, <col>, <val>).

        """

        center_row, center_col, center_val = \
                self.get_left_entry(self.nrows / 2, self.ncols / 2)

        return (center_row, center_col, center_val)

    def get_left_entry(self, row, col):
        """Get the entry to the left of a specified entry, skipping over bulges.

        Arguments:
            row --
            col --

        Returns:
            ...

        """

        new_col = col

        while self.is_stretch(row, new_col):
            new_col -= 1

        return (row, new_col, self.get(row, new_col))

    def get_right_entry(self, row, col):
        """Get the entry to the right of a specified entry, skipping over
        bulges.

        Arguments:
            row --
            col --

        Returns:
            ...

        """

        new_col = col

        while self.is_stretch(row, new_col):
            new_col += 1

        return (row, new_col, self.get(row, new_col))

    def get_neighbouring_entries(self, row, col):
        """Get the neighbouring entries of a specific entry in the beta-matrix.

        Arguments:
            row -- the row index of the entry.
            col -- the column index of the entry.

        Returns:
            a generator over the (<row>, <col>, <val>) entries neighbouring a
            user-specified beta-matrix origin entry.

        """

        # check top neighbour
        if row > 0 and \
           not self.is_padding(row - 1, col) and \
           not self.is_stretch(row - 1, col):
               yield (row - 1, col, self.get(row - 1, col))

        # check bottom neighbour
        if row < self.nrows - 1 and \
           not self.is_padding(row + 1, col) and \
           not self.is_stretch(row + 1, col):
               yield (row + 1, col, self.get(row + 1, col))

        # check left neighbour
        if col > 0:
            if not self.is_padding(row, col - 1):
                yield self.get_left_entry(row, col - 1)

        # check right neighbour
        if col < self.ncols - 1:
            if not self.is_padding(row, col + 1):
                yield self.get_right_entry(row, col + 1)

        return

    def get_alternate_repr(self):
        buf = self.get_header(root_sheet_id=True) + os.linesep
        
        for row in self.matrix:
            buf += "".join(row) + os.linesep

        return buf[:-1]

    def grow(self, start_row=None, start_col=None):
        new_matrix = copy.deepcopy(self)
        new_matrix.clear()

        if start_row == None and start_col == None:
            center_entry = self.get_center_entry()
        else:
            center_entry = self.get(start_row, start_col)

        visited = defaultdict(bool)
        bulges = defaultdict(lambda : None)
        queue = deque()
        queue.append(center_entry)

        while len(queue) > 0:
            row, col, val = queue.popleft()

            if visited[(row, col, val)]:
                continue

            visited[(row, col, val)] = True

            new_matrix.set(row, col, val)

            # add the beta-bulge stretch entries if requried.
            bulge_range = bulges[(row, col, val)]

            if bulge_range != None:
                for j in xrange(bulge_range[0], bulge_range[1]):
                    new_matrix.set(row, j, "-")

            new_matrix.calc_nentries()
            new_matrix.set_sheet_id("%s_s%d"
                                    % (new_matrix.sheet_id, new_matrix.nentries))

            # do not yield matrices with fewer than 3 residues because our
            # smallest feature unit is a trimer which contains 3 residues.
            if new_matrix.nentries > 2:
                yield new_matrix

            for next_entry in self.get_neighbouring_entries(row, col):
                next_row, next_col, next_val = next_entry

                # record the bulge column ranges for addition later.
                if row == next_row and abs(col - next_col) > 1:
                    if col < next_col:
                        bulges[next_entry] = (col + 1, next_col)
                    else:
                        bulges[next_entry] = (next_col + 1, col)

                if not visited[next_entry]:
                    queue.append(next_entry)

        return

    def save_grow_expansions(self, fn=None, folder=None, return_paths=False):
        """Write all of the subqueries "grown" using the self.grow() function
        into a single *.bm file.

        Arguments:
            fn -- path to the file to save each subquery.

        Returns:
            None

        """

        fn_paths = []

        if fn != None:
            f = open(fn, "w")
            for subquery in self.grow():
                f.write(repr(subquery) + os.linesep)
            f.close()
        
        elif folder != None:
            if not os.path.exists(folder):
                os.mkdir(folder)

            for subquery in self.grow():
                fn = os.path.join(folder, "%s.bm" %
                                  subquery.get_sheet_id())
                f = open(fn, "w")
                f.write(repr(subquery) + os.linesep)
                f.close()

                if return_paths:
                    fn_paths.append(fn)

        if return_paths:
            return fn_paths

        return

    def save_class_expansions(self, folder, translations, maxsize=None,
                              return_paths=False, grow=False):
        """Translate a beta-matrix by reducing its alphabet to its residue
        classes.

        Arguments:
            folder -- path to the directory to save each expansion.
            translations --
            maxsize --
            k -- the total number of residues to mutate within their respective
            residue classes (default=None). If None then assume that all
            residues are to be mutated.

        Returns:
            None

        """

        if not os.path.exists(folder):
            os.mkdir(folder)

        if maxsize == None:
            maxsize = self.nentries

        if grow:
            for m in self.grow():
                if m.nentries > maxsize:
                    break

                fn = os.path.join(folder, "%s.bm" % m.get_sheet_id()) 
                m.save_translated_matrices(fn, translations)
        else:
            fn = os.path.join(folder, "%s.bm" % m.get_sheet_id()) 
            m.save_translated_matrices(fn, translations)

        return

    def gen_expansions(self, translations, k=None):
        """Perform an expansion of this beta-matrix, where all k-residue
        combinations are translated to their residue classes.

        Arguments:
            k -- number of residues in the query to be translated
            (default=None). If k is not specified then all the residues in the
            beta-matrix will be translated.

        Returns:
            a generator over nCk query expansions, where n is the total number
            of residues in the query and k is the number of residues in the
            query to be translated.

        """

        matrix_coords = [(row, col) for (row, col, _) in
                         self.gen_entries()]

        if k == None:
            k = self.nentries

        for wildcard_coords in itertools.combinations(matrix_coords, k):
            queue = deque()
            queue.append(self)

            while len(queue) > 0:
                curr_matrix = queue.popleft()
                entry = curr_matrix.is_unexpanded(wildcard_coords)

                if entry == None:
                    curr_matrix.nmutations = k
                    yield curr_matrix
                else:
                    row, col, aa = entry

                    # translate each entry into its residue-class.
                    for res_class in translations.res_to_classes[aa]:
                        new_matrix = copy.deepcopy(curr_matrix)
                        new_matrix.matrix[row][col] = res_class
                        queue.append(new_matrix)

        return

    def save_expansions(self, folder, translations, k=None, overwrite=False):
        """
        Arguments:
            folder -- path to the directory to save each expansion *.bm file.
            translations -- the Translations object to be used.
            k -- the number of tolerated intra-class mutations within a
            beta-matrix (default=None).

        The expansions of each original beta-matrix is saved to a single file
        located in the directory specified by the <folder> argument.

        """

        was_overwritten = False

        if k == None:
            for matrix in self.gen_expansions(translations, k):
                fn = os.path.join(folder,
                                  "%s.bm" % matrix.get_sheet_id())

                if os.path.exists(fn):
                    if overwrite and not was_overwritten:
                        f = open(fn, "w")
                        was_overwritten = True
                    f = open(fn, "a")
                else:
                    f = open(fn, "w")

                f.write(repr(matrix) + os.linesep)
                f.close()
        else:
            for matrix in self.gen_expansions(translations, k):
                fn = os.path.join(folder,
                                  "%s.bm" % matrix.get_sheet_id())
                f = open(fn, "w")

                for enum_matrix in matrix.gen_enumerations(translations):
                    f.write(repr(enum_matrix) + os.linesep)

                f.close()

        return

    def save_translated_matrices(self, fn, translations):
        """
        """

        f = open(fn, "w")

        for trans_matrix in self.gen_translated_matrices(translations):
            f.write(repr(trans_matrix) + os.linesep)

        f.close()

        return

    def gen_entries(self):
        """Create a generator over each non-empty entry (e.g. not '.' nor '-')
        in the beta-matrix. Each entry is a tuple in the form of (<row>, <col>,
        <amino acid>). This method was designed to make it easier to perform
        computations over each amino acid in the beta-matrix.

        Returns:
            a generator over each non-empty (<row>, <col>, <amino acid>)
            entry in this beta-matrix.

        """

        for row in xrange(self.nrows):
            for col in xrange(self.ncols):
                if self.get(row, col) == "." or self.get(row, col) == '-':
                    continue
                yield (row, col, self.get(row, col))

        return

    def is_unexpanded(self, coords):
        """Check if all non-empty entries in a given set of beta-matrix
        coordinates have been translated.

        Returns:
            a (row, col, val) tuple of the first occurrence of a non-empty,
            unexpanded matrix entry from the given set of coordinates;
            otherwise, return None.

        """

        for row, col in coords:
            val = self.get(row, col)
            if val.isupper():
                return row, col, val

        return None

    def next_unenumerated_entry(self):
        """Get the next entry in the beta-matrix that hasn't been enumerated
        from a residue-class character into a single-letter amino acid code yet
        i.e. "unenumerated". Entries are generated left-to-right from the top
        row to the bottom row.

        Arguments:
            None

        Returns:
            a tuple of the form (<row>, <col>, <val>) representing a single
            entry in the beta-matrix if the beta-matrix still contains
            unenumerated entries. Otherwise, a None is returned.

        """

        for row, col, val in self.gen_entries():
            if val.islower():
                return row, col, val

        return None

    def next_amino_acid_entry(self):
        for row, col, val in self.gen_entries():
            if val.isupper():
                return row, col, val

        return None

    def gen_enumerations(self, translations):
        """Translate all the residue class entries in this beta-matrix back into
        their actual amino acids. We refer to this process as "enumeration"
        rather than query expansion. For example:

            [["n", "s"], [".", "C"]]
            
            (this is an example of a "partial expansion")
            
        is the beta-matrix, where the lower-case entries represent residues
        which were translated into their residue classes. The beta-matrix
        enumeration process "reverse translates" this beta-matrix into
        equivalent beta-matrices where each residue class is replaced with its
        containing residues. For example:

            "n" --> { A, G }
            "s" --> { C, P }

        denotes the single-letter codes of the amino acids classed as "n" and
        "s", respectively. Enumerating this beta-matrix will result in the
        following beta-matrices:

            1. [["A", "C"], [".", "C"]]
            2. [["G", "C"], [".", "C"]]
            3. [["A", "P"], [".", "C"]]
            4. [["G", "P"], [".", "C"]]

        I have called this entire process "beta-matrix enumeration", even though
        it is a form of query expansion because we have a couple of different
        types of query expansions in BetaSearch.

        Returns:
            a generator over each enumerated BetaMatrix object.

        """

        queue = deque()
        queue.append(self)

        while len(queue) > 0:
            curr_matrix = queue.popleft()
            entry = curr_matrix.next_unenumerated_entry()
            
            if entry == None:
                yield curr_matrix
            else:
                row, col, res_class = entry

                # reverse translate each residue class into a single-letter
                # amino acid code.
                for aa in translations.classes[res_class]:
                    new_matrix = copy.deepcopy(curr_matrix)
                    new_matrix.set(row, col, aa)
                    queue.append(new_matrix)

        return

    def gen_translated_matrices(self, translations):
        queue = deque()
        queue.append(self)

        while len(queue) > 0:
            curr_matrix = queue.popleft()
            entry = curr_matrix.next_amino_acid_entry()
            
            if entry == None:
                yield curr_matrix
            else:
                row, col, aa = entry

                # reverse translate each residue class into a single-letter
                # amino acid code.
                for res_class in translations.res_to_classes[aa]:
                    new_matrix = copy.deepcopy(curr_matrix)
                    new_matrix.set(row, col, res_class)
                    queue.append(new_matrix)

        return


def beta_matrices_from_file(bm_fn):
    """Create a generator over the beta-matrices in a *.bm file.

    Arguments:
        bm_fn -- path to the *.bm which contains ascii representations of
                 beta-matrices.

    Returns:
        a generator over the beta-matrices in a *.bm file as BetaMatrix objects.

    """

    matrix_lines = []
    i = 0

    f = open(bm_fn, "r")

    for line in f:
        line = line.strip()
        yield BetaMatrix(line=line)

    f.close()

    return


def beta_matrices_from_dir(folder, pattern):
    """Create a generator over the beta-matrices stored in each *.bm in a given
    directory.

    Arguments:
        folder -- path to the directory containing all the *.bm files.
        pattern -- the glob pattern of each *.bm file (the file extension need
                   not even be *.bm).

    Returns:
        None

    """

    for fn in glob.iglob(os.path.join(folder, pattern)):
        yield BetaMatrix(fn)

    return


def trimers_from_file(trimers_fn):
    """Parse the contents of <trimers_fn> into a dictionary where

        dict[t->num] --> t_id
        dict[0] --> "1,LWT"

    this function is intended to be used during the ranking phase after
    verification where we score the trimer similaries between a matching
    beta-matrix and the query.

    Returns:
        a dictionary.

    """

    trimers = {}

    f = open(trimers_fn)

    for line in f:
        line = line.strip()
        for t in line.split():
            t_vals = t.split(":")
            t_id = t_vals[2]
            t_num = int(t_vals[-1])
            trimers[t_num] = t_id

    f.close()

    return trimers


def parse_beta_matrix_header(bm_line):
    """Parse a string representation of a beta-matrix for its header contents,
    storing each field value in a dictionary for ease-of-use.

    """

    vals = bm_line.strip().split(":")
    header = { "sheet_id"      : ":".join(vals[:2]),
               "root_sheet_id" : vals[0],
               "exp_num"       : int(vals[1]),
               "n_rows"        : int(vals[2]),
               "n_cols"        : int(vals[3]),
               "is_barrel"     : vals[4],
               "n_entries"     : int(vals[5]),
               "lines"         : list(list(row.strip()) for row in
                                      vals[-1].split(","))}

    return header


if __name__ == "__main__":
    b = BetaMatrix("test-query.txt")
    t = Translations("res-classes.txt")

    # b.save_enum_expansions("./test-expansions", translations=t, maxsize=4, k=1)
    b.save_grow_expansions(folder="./test-queries")

    """
    b.save_grown_queries("1HHGA-subqueries.txt")

    for m in read_beta_matrices("1HHGA-subqueries.txt"):
        print m
    """

    """
    for m in read_beta_matrices("1HHGA-query.txt"):
        print m
    """

    """
    q = BetaMatrix("test-query.txt")
    t = Translations("res-classes.txt")

    for expanded in q.gen_expansions(t, k=1):
        pprint(expanded.matrix)
        print
    """

    """
    for expanded in q.expansions():
        pprint(expanded.matrix)
        for enumerated in expanded.enumerations(q.translations):
            print "\t", enumerated.matrix
    """

    """
    pprint(q.beta_matrix.matrix)
    row, col, entry = q.beta_matrix.get_center_entry()
    for x in q.beta_matrix.get_neighbouring_entries(row, col):
        print x
    """

    """
    for mat in q.beta_matrix.grow():
        pprint(mat.matrix)
    """
