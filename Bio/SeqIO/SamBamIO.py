# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the SAM/BAM file format.

The SAM (text) and BAM (binary) file format is described here:
http://samtools.sourceforge.net/SAM1.pdf

You are expected to use this module via the Bio.SeqIO functions under
the format names "sam-ref" and "bam-ref" for the possibly annotated
reference sequences in SAM or BAM respectively. SAM v1.5 defines how
annotation is done with dummy reads, and how to embedded the reference
sequence as a dummy read.

>>> from Bio import SeqIO
>>> for ref in SeqIO.parse("SamBam/ex1.sam", "sam-ref"):
...     print ref.id

Or, using the BAM format should behave exactly the same:

>>> from Bio import SeqIO
>>> for ref in SeqIO.parse("SamBam/ex1.bam", "bam-ref"):
...     print ref.id

You can do things like convert from GenBank to SAM,

>>> SeqIO.convert("GenBank/NC_000932.gb", "gb", "NC_000932.sam", "sam-ref")
1
>>> import os
>>> assert 0 == os.system("samtools view -S -b NC_000932.sam | samtools sort - NC_000932")
>>> assert 0 == os.system("samtools index NC_000932.bam")

Internally this module calls Bio.Sequencing.SamBam which offers a more
SAM/BAM specific interface.
"""

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.Interfaces import SequentialSequenceWriter

from Bio.Sequencing import SamBam

def _hack(iterator):
    for read in iterator:
        if read.flag != 516 or read.qname != read.rname:
            continue
        yield SeqRecord(Seq(read.seq, single_letter_alphabet), id=read.qname)

def SamRefIterator(handle):
    """Generator function to iterate over reference sequences in SAM format."""
    return _hack(SamBam.SamIterator(handle))

def BamRefIterator(handle):
    """Generator function to iterate over reference sequences in BAM format."""
    return _hack(SamBam.BamIterator(handle))

class SamRefWriter(SequentialSequenceWriter):
    #Writing the header is a problem for a single pass strategy...
    #Currently assume one record and/or users can sort it out afterwards.
    def write_record(self, record):
        auto_id = 0
        ref = record.id
        cigar = "%i=" % len(record)
        qual = "*" #Get PHRED quality if available
        co = "CO:Z:%s" % record.description.replace("\t", " ")
        parts = [ref, "516", ref, "1", "0", cigar, "*", "0", "0", str(record.seq), qual, co]
        self.handle.write("@SQ\tSN:%s\tLN:%i\n" % (ref, len(record)))
        self.handle.write("\t".join(parts) + "\n")
        for f in record.features:
            def_flag = 768
            if f.id and f.id != "<unknown id>":
                id = f.id
            elif f.sub_features:
                id = "UNNAMED-TEMP-ID-%i" % auto_id
                auto_id += 1
            else:
                id = "*"
            if f.sub_features:
                f_list = f.sub_features
                def_flag += 0x1 #multipart
                def_flag += 0x2 #all properly mapped
                next_ref = ref
            else:
                f_list = [f]
                next_ref = "*"
            #TODO - First/last flag, and next strand flags
            for sf in f_list:
                assert not sf.sub_features
                if sf.location.strand == -1:
                    flag = str(def_flag + 0x10)
                    seq = str(sf.extract(record.seq).reverse_complement())
                else:
                    flag = str(def_flag)
                    seq = str(sf.extract(record.seq))
                start = str(sf.location.start + 1)
                cigar = "%i=" % len(seq)
                co = "CO:Z:%s feature" % f.type
                parts = [id, flag, ref, start, "30", cigar, next_ref, "0", "0", seq, "*", co]
                self.handle.write("\t".join(parts) + "\n")

def _test():
    """Run the module's doctests (PRIVATE).

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..", "..", "..", "Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..", "..", "..", "Tests"))
        doctest.testmod()
        print "Done"
        os.chdir(cur_dir)
        del cur_dir
    elif os.path.isdir(os.path.join("Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        print "Done"
        os.chdir(cur_dir)
        del cur_dir

if __name__ == "__main__":
    _test()
