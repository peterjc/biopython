# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Additional unit tests for Bio.SeqIO.convert(...) function."""
import os
import unittest
from Bio.Seq import UnknownSeq
from Bio import SeqIO
from Bio.SeqIO import QualityIO
from Bio.Alphabet import generic_protein, generic_nucleotide, generic_dna

class IndexDictTests(unittest.TestCase) :
    """Cunning unit test where methods are added at run time."""
    def simple_check(self, filename, format, alphabet) :
        id_list = [rec.id for rec in \
                   SeqIO.parse(open(filename), format, alphabet)]
        rec_dict = SeqIO.indexed_dict(filename, format, alphabet)
        self.assertEqual(set(id_list), set(rec_dict.keys()))
        #This is redundant, I just want to make sure len works:
        self.assertEqual(len(id_list), len(rec_dict))
        for key in id_list :
            self.assert_(key in rec_dict)
            self.assertEqual(key, rec_dict[key].id)
            
tests = [
    ("Ace/contig1.ace", "ace", generic_dna),
    ("Ace/consed_sample.ace", "ace", None),
    ("Ace/seq.cap.ace", "ace", generic_dna),
    ("Quality/example.fastq", "fastq", None),
    ("Quality/example.fastq", "fastq-sanger", generic_dna),
    #Can't yet index line wrapped FASTQ files...
    #("Quality/tricky.fastq", "fastq", generic_nucleotide),
    ("Quality/sanger_faked.fastq", "fastq-sanger", generic_dna),
    ("Quality/solexa_faked.fastq", "fastq-solexa", generic_dna),
    ("Quality/illumina_faked.fastq", "fastq-illumina", generic_dna),
    ("Embl/U87107.embl", "embl", None),
    ("Embl/TRBG361.embl", "embl", None),
    ("GenBank/NC_000932.faa", "fasta", generic_protein),
    ("GenBank/NC_005816.faa", "fasta", generic_protein),
    ("GenBank/NC_005816.ffn", "fasta", generic_dna),
    ("GenBank/NC_005816.fna", "fasta", generic_dna),
    ("GenBank/NC_005816.gb", "gb", None),
    ("GenBank/cor6_6.gb", "genbank", None),
    ("Phd/phd1", "phd", generic_dna),
    ("Phd/phd2", "phd", None),
    ("Phd/phd_solexa", "phd", generic_dna),
    ("Phd/phd_454", "phd", generic_dna),
    ("NBRF/B_nuc.pir", "pir", generic_nucleotide),
    ("NBRF/Cw_prot.pir", "pir", generic_protein),
    ("NBRF/clustalw.pir", "pir", None),
    ("SwissProt/sp001", "swiss", None),
    ("SwissProt/sp010", "swiss", None),
    ("SwissProt/sp016", "swiss", None),
    ("Roche/E3MFGYR02_random_10_reads.sff", "sff", generic_dna),
    ]
for filename, format, alphabet in tests :
    def funct(fn,fmt,alpha) :
        f = lambda x : x.simple_check(fn, fmt, alpha)
        f.__doc__ = "Index %s file %s" % (fmt, fn)
        return f
    setattr(IndexDictTests, "test_%s_%s" \
            % (filename.replace("/","_").replace(".","_"), format),
            funct(filename, format, alphabet))
    del funct

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
