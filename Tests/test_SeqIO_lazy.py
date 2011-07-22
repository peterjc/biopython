import unittest
from StringIO import StringIO

from Bio import SeqIO
from Bio.Alphabet import generic_protein, generic_nucleotide
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio.SeqIO.FastaIO import LazySeqRecordFasta
from Bio.SeqIO.QualityIO import LazySeqRecordFastqSanger, LazySeqRecordFastqSolexa

class TestSimpleRead(unittest.TestCase):
        
    def test_easy_fasta(self):
        """FASTA single entry."""
        filename = "Fasta/centaurea.nu"
        r = SeqIO.read(filename, "fasta", generic_nucleotide)
        
        h = open(filename, "rU")
        h.read()
        p = LazySeqRecordFasta(h, 0, h.tell(), len(r), slice(None, None, None),
                               r.seq.alphabet)
        self.compare(p, r)

        for s in [slice(5,10),
                  slice(10,20,1),
                  slice(5,-5),
                  ]:
            p1 = p[s]
            r1 = p[s]
            self.compare(p1, r1)
        h.close()

    def test_first_fasta(self):
        """FASTA first entry."""
        filename = "Quality/example.fasta"
        r = SeqIO.parse(filename, "fasta", generic_nucleotide).next()
        
        h = open(filename)
        h.seek(6) #random
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
        h.close()

    def test_first_fastq(self):
        """FASTQ Sanger first entry."""
        filename = "Quality/example.fastq"
        r = SeqIO.parse(filename, "fastq", generic_nucleotide).next()
        
        h = open(filename)
        p = LazySeqRecordFastqSanger(h, 0, None, len(r), slice(None, None, None),
                                     r.seq.alphabet)
        self.compare(p, r)

        for s in [slice(5,10),
                  slice(10,20,1),
                  slice(5,-5),
                  ]:
            p1 = p[s]
            r1 = p[s]
            self.compare(p1, r1)
        self.assertEqual(p.format("fastq"), r.format("fastq"))
        h.close()

    def test_first_fastq_solexa(self):
        """FASTQ Solexa first entry."""
        filename = "Quality/solexa_example.fastq"
        r = SeqIO.parse(filename, "fastq-solexa", generic_nucleotide).next()
        
        h = open(filename)
        h.seek(19) #random
        p = LazySeqRecordFastqSolexa(h, 0, 120, len(r), slice(None, None, None),
                                     r.seq.alphabet)
        self.compare(p, r)
        
        for s in [slice(5,10),
                  slice(10,20,1),
                  slice(5,-5),
                  ]:
            p1 = p[s]
            r1 = p[s]
            self.compare(p1, r1)
        self.assertEqual(p.format("fastq-solexa"), '@SLXA-B3_649_FC8437_R1_1_1_610_79\nGATGTGCAATACCTTTGTAGAGGAA\n+SLXA-B3_649_FC8437_R1_1_1_610_79\nYYYYYYYYYYYYYYYYYYWYWYYSU\n')
        self.assertEqual(r.format("fastq-solexa"), '@SLXA-B3_649_FC8437_R1_1_1_610_79\nGATGTGCAATACCTTTGTAGAGGAA\n+\nYYYYYYYYYYYYYYYYYYWYWYYSU\n')
        h.close()

    def compare(self, p, r):
        self.assertEqual(len(p), len(r))
        self.assertEqual(p.id, r.id)
        self.assertEqual(p.name, r.name)
        self.assertEqual(p.description, r.description)
        self.assertEqual(p.annotations, r.annotations)
        self.assertEqual(p.letter_annotations, r.letter_annotations)
        self.assertEqual(len(p.features), len(r.features))
        self.assertEqual(p.dbxrefs, r.dbxrefs)
        self.assertEqual(str(p.seq), str(r.seq))
        #self.assertEqual(p.format("fasta"), r.format("fasta"))

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
