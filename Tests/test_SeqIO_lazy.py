import unittest
from StringIO import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import LazySeqRecordFasta
from Bio.Alphabet import generic_protein, generic_nucleotide

class TestSimpleRead(unittest.TestCase):
        
    def test_easy(self):
        r = SeqIO.read("Fasta/centaurea.nu", "fasta", generic_nucleotide)
        
        h = open("Fasta/centaurea.nu")
        p = LazySeqRecordFasta(h, 0, None, len(r), slice(None, None, None),
                               r.seq.alphabet)
        self.compare(p, r)

        for s in [slice(5,10),
                  slice(10,20,1),
                  slice(5,-5),
                  ]:
            p1 = p[s]
            r1 = p[s]
            self.compare(p1, r1)
    
    def compare(self, p, r):
        self.assertEqual(len(p), len(r))
        self.assertEqual(p.id, r.id)
        self.assertEqual(p.name, r.name)
        self.assertEqual(p.description, r.description)
        self.assertEqual(p.annotations, r.annotations)
        self.assertEqual(str(p.seq), str(r.seq))

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
