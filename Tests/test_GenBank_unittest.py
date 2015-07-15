# Copyright 2013 by Kai Blin.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest
import warnings
from os import path

from Bio import BiopythonParserWarning
from Bio import GenBank
from Bio import SeqIO


class GenBankTests(unittest.TestCase):
    def test_invalid_product_line_raises_value_error(self):
        """Test GenBank parsing invalid product line raises ValueError"""
        def parse_invalid_product_line():
            rec = SeqIO.read(path.join('GenBank', 'invalid_product.gb'),
                             'genbank')
        self.assertRaises(ValueError, parse_invalid_product_line)

    def test_genbank_read(self):
        with open(path.join("GenBank", "NC_000932.gb")) as handle:
            record = GenBank.read(handle)
        self.assertEqual(['NC_000932'], record.accession)

    def test_genbank_read_multirecord(self):
        with open(path.join("GenBank", "cor6_6.gb")) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    def test_genbank_read_invalid(self):
        with open(path.join("GenBank", "NC_000932.faa")) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    def test_genbank_read_no_origin_no_end(self):
        with open(path.join("GenBank", "no_origin_no_end.gb")) as handle:
            self.assertRaises(ValueError, GenBank.read, handle)

    # Evil hack with 000 to manipulate sort order to ensure this is tested
    # first (otherwise something silences the warning)
    def test_000_genbank_bad_loc_wrap_warning(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            with open(path.join("GenBank", "bad_loc_wrap.gb")) as handle:
                # self.assertRaises(BiopythonParserWarning, GenBank.read, handle)
                try:
                    record = GenBank.read(handle)
                except BiopythonParserWarning as e:
                    self.assertEqual(str(e), "Non-standard feature line wrapping (didn't break on comma)?")
                else:
                    self.assertTrue(False, "Expected specified BiopythonParserWarning here.")

    def test_genbank_bad_loc_wrap_parsing(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonParserWarning)
            with open(path.join("GenBank", "bad_loc_wrap.gb")) as handle:
                record = GenBank.read(handle)
                self.assertEqual(1, len(record.features))
                loc = record.features[0].location
                self.assertEqual(loc, "join(3462..3615,3698..3978,4077..4307,4408..4797,4876..5028,5141..5332)")

    def test_multi_exon_rev_strand(self):
        record = SeqIO.read("GenBank/NC_000932_fragment.gb", "gb")
        # Two versions of same CDS, using join(complement(...),complement(...))
        # versus complement(join(...,...))
        self.assertEqual(record.features[1].qualifiers["gene"], ["rps16"])
        self.assertEqual(record.features[2].qualifiers["gene"], ["rps16_alt"])
        self.assertEqual(str(record.features[1].location), str(record.features[2].location))
        for feature in record.features:
            if feature.type == "CDS":
                gene = feature.qualifiers['gene'][0]
                expected = feature.qualifiers['translation'][0]
                table = feature.qualifiers['transl_table'][0]
                trans = str(feature.extract(record.seq).translate(table, cds=True))
                self.assertTrue(trans == expected or trans == expected + "*",
                                "Translation problem for %s, have %r not %r"
                                % (gene, trans, expected))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
