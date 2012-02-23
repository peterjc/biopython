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
try:
    from itertools import izip_longest
except ImportError:
    #That was renamed in Python 3, should be fixed by 2to3 but not:
    #http://bugs.python.org/issue11438
    try:
        from itertools import zip_longest as izip_longest
    except ImportError:
        #Not present at all on Python 2.5 or Jython
        def izip_longest(A, B):
            a_iter = iter(A)
            b_iter = iter(B)
            while True:
                try:
                    a = a_iter.next()
                except StopIteration:
                    a = None
                try:
                    b = b_iter.next()
                except StopIteration:
                    b = None
                if a is None and b is None:
                    break
                else:
                    yield a, b

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
        #Note considerable overlap with method in test_SamBam_pysam.py
        self.assertEqual(a_iter.nreferences, b_iter.nreferences)
        self.assertEqual(a_iter.references, b_iter.references)
        self.assertEqual(a_iter.lengths, b_iter.lengths)
        if a_iter.header and b_iter.header:
            self.assertEqual(a_iter.header, b_iter.header)
        for a, b in izip_longest(a_iter, b_iter):
            self.assertFalse(b is None, "Extra read in a:\n%s" % a)
            self.assertFalse(a is None, "Extra read in b:\n%s" % b)
            self.assertEqual(a.qname, b.qname)
            self.assertEqual(a.flag, b.flag,
                             "%r vs %r for:\n%s\n%s"  % (a.flag, b.flag, a, b))
            self.assertEqual(a.rname, b.rname,
                             "%r vs %r for:\n%s\n%s"  % (a.rname, b.rname, a, b))
            self.assertEqual(a.pos, b.pos,
                             "%r vs %r for:\n%s\n%s"  % (a.pos, b.pos, a, b))
            self.assertEqual(a.alen, b.alen,
                             "%r vs %r for:\n%s\n%s"  % (a.alen, b.alen, a, b))
            self.assertEqual(a.aend, b.aend,
                             "%r vs %r for:\n%s\n%s"  % (a.aend, b.aend, a, b))
            self.assertEqual(a.mapq, b.mapq)
            self.assertEqual(a.cigar, b.cigar,
                             "%r vs %r for:\n%s\n%s" % (a.cigar, b.cigar, a, b))
            self.assertEqual(a.rnext, b.rnext,
                             "%r vs %r for:\n%s\n%s" % (a.rnext, b.rnext, a, b))
            self.assertEqual(a.mpos, b.mpos)
            self.assertEqual(a.isize, b.isize)
            self.assertEqual(a.seq, b.seq,
                             "%r vs %r for:\n%s\n%s" % (a.seq, b.seq, a, b))
            self.assertEqual(a.qual, b.qual)
            self.assertEqual(a.tags.keys(), b.tags.keys())
            for key in a.tags:
                self.assertEqual(a.tags[key][0], b.tags[key][0],
                       "%s tag %s has codes %s vs %s" % \
                       (a.qname, key, a.tags[key][0], b.tags[key][0]))
                if a.tags[key][0]=="f":
                    self.assertTrue(str(a.tags[key][1]) == str(b.tags[key][1]) or \
                                    abs(a.tags[key][1] - b.tags[key][1]) < 0.00001, \
                                    "%s tag %s have %f vs %f" % \
                                    (a.qname, key, a.tags[key][1], b.tags[key][1]))
                elif a.tags[key][0]=="Bf":
                    pass
                else:
                    self.assertEqual(a.tags[key][1], b.tags[key][1],
                           "%s tag %s:%s has %s vs %s" % \
                           (a.qname, key, a.tags[key][0], a.tags[key][1], b.tags[key][1]))
            #Float formating in tags is annoying...
            #assert str(a) == str(b), "Reads disagree,\n%s\n%s\n" % (a,b)

    def test_bam_vs_bam_ex1(self):
        self.compare(BamIterator(open("SamBam/ex1_header.bam", "rb")),
                     BamIterator(gzip.open("SamBam/ex1_header.bam"), gzipped=False))

    def test_sam_vs_bam_ex1_header(self):
        self.compare(SamIterator(open("SamBam/ex1_header.sam")),
                     BamIterator(open("SamBam/ex1_header.bam", "rb")))
        self.compare(SamIterator(open("SamBam/ex1_header.sam")),
                     BamIterator(gzip.open("SamBam/ex1_header.bam"), gzipped=False))

    def test_sam_vs_bam_bins(self):
        self.compare(SamIterator(open("SamBam/bins.sam")),
                     BamIterator(open("SamBam/bins.bam", "rb")))
        self.compare(SamIterator(open("SamBam/bins.sam")),
                     BamIterator(gzip.open("SamBam/bins.bam"), gzipped=False))

    def test_sam_vs_bam_tags(self):
        self.compare(SamIterator(open("SamBam/tags.sam")),
                     BamIterator(open("SamBam/tags.bam", "rb")))
        self.compare(SamIterator(open("SamBam/tags.sam")),
                     BamIterator(gzip.open("SamBam/tags.bam"), gzipped=False))

    def test_tags_sam(self):
        for read in SamIterator(open("SamBam/tags.sam")):
            #TODO - Test API for getting tag values
            tag = str(read).rstrip("\n").split("\t")[-1]
            self.assertEqual(read.qname, "tag_" + tag,
                             "%s vs tag of %s" % (read.qname, tag))

    def test_tags_bam(self):
        for read in BamIterator(open("SamBam/tags.bam", "rb")):
            tag = str(read).rstrip("\n").split("\t")[-1]
            self.assertTrue(read.qname.startswith("tag_" + tag[:5]))
            if read.qname == "tag_" + tag:
                continue
            if ":f:" in tag:
                old = float(read.qname.split(":")[2])
                new = float(tag.split(":")[2])
                self.assertTrue(_comp_float(old, new),
                                "%s vs tag of %s" % (read.qname, tag))
            elif ":B:" in tag:
                self.assertTrue(read.qname.startswith("tag_" + tag[:5]),
                                "%s vs tag of %s" % (read.qname, tag))
                if ":B:f," in read.qname:
                    old = tuple(float(v) for v in read.qname.split(":B:f,")[1].split(","))
                    new = tuple(float(v) for v in tag.split(":")[2][2:].split(","))
                    self.assertEqual(len(old), len(new),
                                     "Ccount mismatch %i vs %i in %s\n%s\n%s\n" \
                                     % (len(old), len(new), read.qname, old, new))
                    for a,b in zip(old, new):
                        self.assertTrue(_comp_float(a, b),
                                        "Mismatch in %s,\n%s\n%s" % (read.qname, old, new))
                else:
                    #Integers
                    old = [int(v) for v in read.qname.split(":B:")[1][2:].split(",")]
                    new = [int(v) for v in tag.split(":")[2][2:].split(",")]
                    self.assertEqual(old, new,
                            "Mismatch in read %s vs %s\n%r\n%r" % (read.qname, tag, old, new))
            else:
                self.assertEqual(read.qname, "tag_" + tag,
                                 "Please check %s vs tag of %s" % (read.qname, tag))


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
        #Idea behind using batches of 16 is for easier cross referencing
        #with a tools like hexdump -C
        for i in range(0, len(old), 16):
            self.assertEqual(old[i:i+16], new[i:i+16],
                             "Mismatch recreating %s bytes %i to %i (0x%x to 0x%x),\n%r\n%r\n" \
                             % (name, i, i+16, i, i+16, old[i:i+16], new[i:i+16]))
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
