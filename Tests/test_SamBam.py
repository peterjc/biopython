# Copyright 2010-2012 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Cross checking our SAM/BAM support against pysam."""

import unittest
import gzip
import os
import sys
from itertools import izip_longest

try:
    #This is in Python 2.6+, but we need it on Python 3
    from io import BytesIO
except ImportError:
    #Must be on Python 2.5 or older
    from StringIO import StringIO as BytesIO

from Bio.Sequencing.SamBam import SamIterator, BamIterator
from Bio.Sequencing.SamBam import SamWriter, BamWriter

def _comp_float(a,b):
    return repr(a)==repr(b) or a==b or abs(float(a)-float(b))<0.000001,

class MiscTests(unittest.TestCase):
    """Misc SAM/BAM tests."""

    def compare(self, a_iter, b_iter):
        assert a_iter.nreferences == b_iter.nreferences, \
            "%r vs %r" % (a_iter.nreferences, b_iter.nreferences)
        assert a_iter.references == b_iter.references, \
            "%r vs %r" % (a_iter.references, b_iter.references)
        assert a_iter.lengths == b_iter.lengths, \
            "%r vs %r" % (a_iter.lengths, b_iter.lengths)
        if a_iter.header and b_iter.header:
            assert a_iter.header == b_iter.header, \
                "%r vs %r" % (a_iter.header, b_iter.header)
        for a, b in izip_longest(a_iter, b_iter):
            assert b is not None, "Extra read in a: %r" % str(a)
            assert a is not None, "Extra read in b: %r" % str(b)
            assert a.qname == b.qname, "%r vs %r" % (a.qname, b.qname)
            assert a.flag == b.flag, "%r vs %r" % (a.flag, b.flag)
            assert a.rname == b.rname, "%r vs %r" % (a.rname, b.rname)
            assert a.pos == b.pos, "%r vs %r" % (a.pos, b.pos)
            assert a.mapq == b.mapq, "%r vs %r" % (a.mapq, b.mapq)
            assert a.cigar == b.cigar, "%r vs %r" % (a.cigar, b.cigar)
            assert a.mrnm == b.mrnm, "%r vs %r" % (a.mrnm, b.mrnm)                                                        
            assert a.mpos == b.mpos, "%r vs %r" % (a.pos, b.pos)
            assert a.isize == b.isize, "%r vs %r" % (a.isize, b.isize)
            assert a.seq == b.seq, "%r vs %r" % (a.seq, b.seq)
            assert a.qual == b.qual, "%r vs %r" % (a.qual, b.qual)
            assert a.tags.keys() == b.tags.keys(), "%r vs %r" % (a.tags.keys(), b.tags.keys())
            for key in a.tags:
                assert a.tags[key][0] == b.tags[key][0], \
                       "%s tag %s has codes %s vs %s" % \
                       (a.qname, key, a.tags[key][0], b.tags[key][0])
                if a.tags[key][0]=="f":
                    assert str(a.tags[key][1]) == str(b.tags[key][1]) or \
                           abs(a.tags[key][1] - b.tags[key][1]) < 0.00001, \
                           "%s tag %s have %f vs %f" % \
                           (a.qname, key, a.tags[key][1], b.tags[key][1])
                elif a.tags[key][0]=="Bf":
                    pass
                else:
                    assert a.tags[key][1] == b.tags[key][1], \
                           "%s tag %s:%s has %s vs %s" % \
                           (a.qname, key, a.tags[key][0], a.tags[key][1], b.tags[key][1])
            #Float formating in tags is annoying...
            #assert str(a) == str(b), "Reads disagree,\n%s\n%s\n" % (a,b)

    def test_reproduce_tags(self):
        self.compare(SamIterator(open("SamBam/tags.sam")),
                     BamIterator(open("SamBam/tags.bam", "rb")))
        self.compare(SamIterator(open("SamBam/tags.sam")),
                     BamIterator(gzip.open("SamBam/tags.bam"), gzipped=False))

    def test_tags_sam(self):
        for read in SamIterator(open("SamBam/tags.sam")):
            #TODO - Test API for getting tag values
            tag = str(read).rstrip("\n").split("\t")[-1]
            assert read.qname == "tag_" + tag, \
                   "%s vs tag of %s" % (read.qname, tag)

    def test_tags_bam(self):
        for read in BamIterator(open("SamBam/tags.bam", "rb")):
            tag = str(read).rstrip("\n").split("\t")[-1]
            assert read.qname.startswith("tag_" + tag[:5])
            if read.qname == "tag_" + tag:
                continue
            if ":f:" in tag:
                old = float(read.qname.split(":")[2])
                new = float(tag.split(":")[2])
                assert _comp_float(old, new), \
                       "%s vs tag of %s" % (read.qname, tag)
            elif ":B:" in tag:
                assert read.qname.startswith("tag_" + tag[:5]), \
                       "%s vs tag of %s" % (read.qname, tag)
                if ":B:f," in read.qname:
                    old = tuple(float(v) for v in read.qname.split(":B:f,")[1].split(","))
                    new = tuple(float(v) for v in tag.split(":")[2][2:].split(","))
                    assert len(old) == len(new), \
                           "Count mismatch %i vs %i in %s\n%s\n%s\n" \
                           % (read.qname, len(old), len(new), old, new)
                    for a,b in zip(old, new):
                        assert _comp_float(a, b), \
                               "Mismatch in %s,\n%s\n%s" % (read.qname, old, new)
                else:
                    #Integers
                    old = [int(v) for v in read.qname.split(":B:")[1][2:].split(",")]
                    new = [int(v) for v in tag.split(":")[2][2:].split(",")]
                    assert old == new, "Mismatch in read %s vs %s\n%r\n%r" % (read.qname, tag, old, new)
            else:
                assert read.qname == "tag_" + tag, \
                       "Please check %s vs tag of %s" % (read.qname, tag)


class RecreateTests(unittest.TestCase):
    """Check can recreate sample files.

    For BAM files, the BGZF block break points are arbitrary. Older
    versions of samtools write full blocks, newer versions try to avoid
    splitting reads between blocks. We therefore compare the 'naked'
    uncompressed binary data from BAM files.
    """

    def check_old_new(self, old, new, name):
        self.assertEqual(old[:100], new[:100])
        self.assertEqual(len(old), len(new))
        #Idea behind using batches of 32 is for easier cross referencing
        #with a tools like hexdump
        for i in range(0, len(old), 32):
            self.assertEqual(old[i:i+32], new[i:i+32],
                             "Mismatch recreating %s bytes %i to %i,\n%r\n%r\n" \
                             % (name, i, i+32, old[i:i+32], new[i:i+32]))
        self.assertEqual(old, new, "Could not recreate %s perfectly" % name)

    def check_bam_to_bam(self, bam_filename):
        handle = BytesIO()
        count = BamWriter(handle, BamIterator(open(bam_filename, "rb")),
                          gzipped=False)
        self.check_old_new(gzip.open(bam_filename).read(),
                           handle.getvalue(), bam_filename)

    def check_sam_to_bam(self, sam_filename, bam_filename):
        handle = BytesIO()
        count = BamWriter(handle, SamIterator(open(sam_filename)),
                          gzipped=False)
        self.check_old_new(gzip.open(bam_filename).read(),
                           handle.getvalue(), bam_filename)

    #TODO - check SAM to SAM (floating points make it tricky)

    def test_tags_b2b(self):
        self.check_bam_to_bam("SamBam/tags.bam")

    def test_tags_s2b(self):
        self.check_sam_to_bam("SamBam/tags.sam", "SamBam/tags.bam")

    def test_ex1_header_b2b(self):
        self.check_bam_to_bam("SamBam/ex1_header.bam")

    def test_ex1_header_s2b(self):
        self.check_sam_to_bam("SamBam/ex1_header.sam", "SamBam/ex1_header.bam")

    def test_bins_b2b(self):
        self.check_bam_to_bam("SamBam/bins.bam")

    def test_bins_s2b(self):
        self.check_sam_to_bam("SamBam/bins.sam", "SamBam/bins.bam")

if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
