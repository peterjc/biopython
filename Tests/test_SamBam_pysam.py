# Copyright 2010-2012 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Cross checking our SAM/BAM support against pysam."""

try:
    import pysam
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install pysam if you want to test against it.")

import unittest
from itertools import izip_longest

from Bio.Sequencing.SamBam import SamIterator, BamIterator

class CrossCheckParsing(unittest.TestCase):
    """Cross checking our SAM/BAM against pysam."""

    def compare(self, a_iter, b_iter):
        #Note this differs from the comparison in test_SamBam.py
        #in order to account for differences with pysam
        #Also it copes with lists (useful with fetch testing)
        if isinstance(a_iter, list) or isinstance(b_iter, list):
            self.assertEqual(len(a_iter), len(b_iter))
        elif isinstance(a_iter, list) or isinstance(b_iter, list):
            pass
        else:
            self.assertEqual(a_iter.nreferences, b_iter.nreferences)
            self.assertEqual(a_iter.references, b_iter.references)
            self.assertEqual(a_iter.lengths, b_iter.lengths)
            if a_iter.header and b_iter.header:
                #pysam doesn't infer a minimal SAM header from BAM header
                self.assertEqual(a_iter.header, b_iter.header)
        for a, b in izip_longest(a_iter, b_iter):
            self.assertFalse(b is None, "Extra read in a:\n%s" % a)
            self.assertFalse(a is None, "Extra read in b:\n%s" % b)
            self.assertEqual(a.qname, b.qname)
            self.assertEqual(a.flag, b.flag,
                             "%r vs %r for:\n%s\n%s"  % (a.flag, b.flag, a, b))
            #See http://code.google.com/p/pysam/issues/detail?id=25
            #self.assertEqual(a.rname, b.rname)
            self.assertEqual(a.pos, b.pos,
                             "%r vs %r for:\n%s\n%s"  % (a.pos, b.pos, a, b))
            #pysam 0.6 uses samtools 0.1.18 where bam_calend doesn't support CIGAR X/=
            #self.assertEqual(a.alen, b.alen,
            #                 "%r vs %r for:\n%s\n%s"  % (a.alen, b.alen, a, b))
            #self.assertEqual(a.aend, b.aend,
            #                 "%r vs %r for:\n%s\n%s"  % (a.aend, b.aend, a, b))
            self.assertEqual(a.mapq, b.mapq)
            self.assertEqual(a.cigar, b.cigar,
                             "%r vs %r for:\n%s\n%s" % (a.cigar, b.cigar, a, b))
            #See http://code.google.com/p/pysam/issues/detail?id=25
            #self.assertEqual(a.mrnm, b.mrnm)
            self.assertEqual(a.mpos, b.mpos)
            self.assertEqual(a.isize, b.isize) #legacy alias
            self.assertEqual(a.seq, b.seq,
                             "%r vs %r for:\n%s\n%s" % (a.seq, b.seq, a, b))
            self.assertEqual(a.qual, b.qual)
            #TODO - tags (not represented same way at the moment)

            #assert str(a) == str(b), "Reads disagree,\n%s\n%s\n" % (a,b)
            #Would compare str(a)==str(b) but pysam does not attempt
            #to return a valid SAM record like this.
            #See http://code.google.com/p/pysam/issues/detail?id=74
            #and http://code.google.com/p/pysam/issues/detail?id=75
    
    #TODO, use pysam on the SAM file
    #http://code.google.com/p/pysam/issues/detail?id=73
    #compare(SamIterator(open("SamBam/ex1.sam")),
    #        pysam.Samfile("SamBam/ex1.sam", "r"))
    #Avoid this by using a copy of the SAM file with a header:

    def test_ex1_header_sam_vs_sam(self):
        self.compare(SamIterator(open("SamBam/ex1_header.sam")),
                     pysam.Samfile("SamBam/ex1_header.sam", "r"))

    def test_ex1_header_sam_vs_bam(self):
        self.compare(SamIterator(open("SamBam/ex1_header.sam")),
                     pysam.Samfile("SamBam/ex1_header.bam", "rb"))

    def test_ex1_header_bam_vs_bam(self):
        self.compare(BamIterator(open("SamBam/ex1_header.bam", "rb")),
                     pysam.Samfile("SamBam/ex1_header.bam", "rb"))

    #pysam:
    #ValueError: file header is empty (mode='r') - is it SAM/BAM format?
    #def test_ex1_bam_vs_bam(self):
    #    self.compare(BamIterator(open("SamBam/ex1.bam", "rb")),
    #                 pysam.Samfile("SamBam/ex1.bam", "rb"))
    #
    #def test_ex1_sam_vs_sam(self):
    #    self.compare(SamIterator(open("SamBam/ex1.sam")),
    #                 pysam.Samfile("SamBam/ex1.sam", "r"))

    def test_tags_sam_vs_sam(self):
        self.compare(SamIterator(open("SamBam/tags.sam")),
                     pysam.Samfile("SamBam/tags.sam", "r"))

    def test_tags_bam_vs_sam(self):
        self.compare(BamIterator(open("SamBam/tags.bam", "rb")),
                     pysam.Samfile("SamBam/tags.sam", "r"))

    def test_tags_sam_vs_bam(self):
        self.compare(SamIterator(open("SamBam/tags.sam")),
                     pysam.Samfile("SamBam/tags.bam", "rb"))

    def test_tags_bam_vs_bam(self):
        self.compare(BamIterator(open("SamBam/tags.bam", "rb")),
                     pysam.Samfile("SamBam/tags.bam", "rb"))


    def test_bins_sam_vs_sam(self):
        self.compare(SamIterator(open("SamBam/bins.sam")),
                     pysam.Samfile("SamBam/bins.sam", "r"))

    def test_bins_bam_vs_sam(self):
        self.compare(BamIterator(open("SamBam/bins.bam", "rb")),
                     pysam.Samfile("SamBam/bins.sam", "r"))

    def test_bins_sam_vs_bam(self):
        self.compare(SamIterator(open("SamBam/bins.sam")),
                     pysam.Samfile("SamBam/bins.bam", "rb"))

    def test_bins_bam_vs_bam(self):
        self.compare(BamIterator(open("SamBam/bins.bam", "rb")),
                     pysam.Samfile("SamBam/bins.bam", "rb"))


    def index_check(self, bam, bai, regions):
        #print bam
        pysam_bam = pysam.Samfile(bam, "rb")
        handle = open(bam, "rb")
        biopy_bam = BamIterator(handle, bai_filename = bai)
        for ref, start, end in regions:
            #print "%s region %i:%i" % (ref, start, end)
            pysam_list = list(pysam_bam.fetch(ref, start, end))
            #print "%s region %i:%i has %i reads (pysam)" \
            #      % (ref, start, end, len(pysam_list))
            biopy_list = list(biopy_bam.fetch(ref, start, end))
            #print "%s region %i:%i has %i reads (biopy)" \
            #      % (ref, start, end, len(biopy_list))
            self.compare(pysam_list, biopy_list)
        handle.close()

    def test_index_ex1(self):
        self.index_check("SamBam/ex1.bam", "SamBam/ex1.bam.bai", [
                            ("chr1", 120, 130),
                        ])

    def test_index_bins(self):
        self.index_check("SamBam/bins.bam", "SamBam/bins.bam.bai", [
                            ("large", 120, 130),
                            ("large", 1000, 1010),
                        ])


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
