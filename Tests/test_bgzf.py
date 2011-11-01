# Copyright 2010-2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test code for working with BGZF files (used in BAM files).

See also the doctests in bgzf.py which are called via run_tests.py
"""

import unittest
import gzip

from Bio.Sequencing.SamBam.bgzf import BgzfWriter, BgzfBlocks

class BgzfTests(unittest.TestCase):
    def test_one(self):
        """Reproduce BGZF compression for BAM file"""
        temp_file = "temp.bgzf"

        #Note this example is from an old version of samtools
        #and all the blocks are full (except the last one)
        h = open("SamBam/ex1.uncompressed.bam", "rb")
        data = h.read()
        h.close()

        h = BgzfWriter(temp_file, "wb")
        h.write(data)
        h.flush()
        h.flush() #Force a dummy empty BGZF block as EOF marker
        h.close()

        h = gzip.open(temp_file)
        new_data = h.read()
        h.close()

        self.assert_(new_data, "Empty BGZF file?")
        self.assertEqual(len(data), len(new_data))
        self.assertEqual(data, new_data)

        h = open("SamBam/ex1.bam", "rb")
        old = list(BgzfBlocks(h))
        h.close()

        h = open(temp_file, "rb")
        new = list(BgzfBlocks(h))
        h.close()

        self.assertEqual(old, new)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
