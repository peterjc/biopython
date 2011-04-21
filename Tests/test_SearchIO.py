# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest

from Bio import SearchIO

class CheckNoExceptions(unittest.TestCase):

    def test_blast_text(self):
        """Parsing BLAST pairwise text"""
        for f in os.listdir("Blast"):
            if f.startswith("bt") and f.endswith(".txt"):
                filename = os.path.join("Blast", f)
                data = open(filename).read()
                if "<HTML>" in data or "<pre>" in data:
                    continue
                try:
                    for query_result in SearchIO.parse(filename, "blast-text"):
                        pass
                except ValueError, e:
                    if "This consumer doesn't handle PSI-BLAST data" in str(e):
                        pass
                    else:
                        raise ValueError("%s - %s" % (f, e))

    def test_blast_xml(self):
        """Parsing BLAST XML"""
        for f in os.listdir("Blast"):
            if f.endswith(".xml"):
                filename = os.path.join("Blast", f)
                for query_result in SearchIO.parse(filename, "blast-xml"):
                    pass

    
    def test_blast_stdtab(self):
        """Parsing BLAST tabular"""
        for f in os.listdir("Blast"):
            if f.endswith(".tabular"):
                filename = os.path.join("Blast", f)
                for query_result in SearchIO.parse(filename, "blast-stdtab"):
                    pass

    def test_fasta_m10(self):
        """Parsing FASTA -m 10"""
        for f in os.listdir("Fasta"):
            if f.endswith(".m10"):
                filename = os.path.join("Fasta", f)
                for query_result in SearchIO.parse(filename, "fasta-m10"):
                    pass

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
