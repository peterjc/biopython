# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys
import subprocess
import unittest

from Bio import MissingExternalDependencyError

from Bio.Sequencing.SamBam.bai import idxstats

#################################################################

#Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

samtools_exe = None
if sys.platform=="win32":
    raise MissingExternalDependencyError("This test does not yet run on Windows")
else:
    import commands
    output = commands.getoutput("samtools")
    #Since "not found" may be in another language, try and be sure this is
    #really the MUSCLE tool's output
    if "not found" not in output \
    and "Tools for alignments in the SAM format" in output:
        samtools_exe = "samtools"

if not samtools_exe:
    raise MissingExternalDependencyError(\
        "Install samtools if you want to run this test.")

#################################################################

class SamtoolsIdxstats(unittest.TestCase):
    def compare(self, bam_filename):
        bai_filename = bam_filename + ".bai"
        #TODO - Switch to subprocess for testing under Windows
        old = commands.getoutput("%s idxstats %s" % (samtools_exe, bam_filename)).strip().split("\n")
        new = ["%s\t%s\t%i\t%i" % values for values in idxstats(bam_filename, bai_filename)]
        self.assertEqual(len(old), len(new))
        for o,n in zip(old, new):
            self.assertEqual(o, n,
                "samtools idxstats %s:\n%s\n\nVersus:\n%s" % (bam_filename, old, new))

    def test_ex1(self):
        self.compare("SamBam/ex1.bam")

    def test_ex1_header(self):
        self.compare("SamBam/ex1_header.bam")

    def test_bins(self):
        self.compare("SamBam/bins.bam")

    def test_tags(self):
        self.compare("SamBam/tags.bam")


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
