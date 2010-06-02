# Copyright 2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Additional unit tests for Bio.SeqIO.SamBamIO (covering SAM and BAM)."""
import os
import unittest
import warnings
from Bio.Alphabet import generic_dna
from Bio.SeqIO import SamBamIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO
from Bio.Data.IUPACData import ambiguous_dna_letters, ambiguous_rna_letters

from seq_tests_common import compare_records

class TestSamVsBam(unittest.TestCase):
    def test_agree(self):
        """Test SAM and BAM parsers agree"""
        #Parse the whole file, but only compare the first 20 to make this fast
        sam = list(SeqIO.parse("SamBam/ex1.sam", "sam"))[:20]
        bam = list(SeqIO.parse("SamBam/ex1.bam", "bam"))[:20]
        self.assert_(compare_records(sam, bam))

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
