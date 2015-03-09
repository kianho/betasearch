###############################################################################
#
# genseqid.py - simple class to generate sequential id numbers
#
# File:    genseqid.py
# Author:  Alex Stivala
# Created: February 2008
#
# $Id: genseqid.py,v 5bbae4254475 2012/05/05 06:48:47 hohkhkh1 $
#
# This class is used to generate sequential id numbers.
#
###############################################################################

class GenSeqId:
    """
    GenSeqId is a class to generate sequential id numbers.
    Just initialize with the first number to use, then call next() to
    get next number in sequence.
    """

    def __init__(self, firstnum = 0):
        """
        Initialize the sequential id generator.

        Parameters:
            firstnum - (default 0), integer to start at
        """
        self.nextnum = firstnum

    def next(self):
        """
        Return the next number in the sequence.
        """
        num = self.nextnum
        self.nextnum += 1
        return num

    
