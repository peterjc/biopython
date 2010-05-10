# Copyright 2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This module is for reading SAM (Sequence Alignment/Map) format files, and
# the closely related binary variant BAM format files, as SeqRecord objects.
# This module is expected to be used via the Bio.SeqIO API.

"""Bio.SeqIO support for the SAM and BAM file formats.

Note that you are expected to use this code via the Bio.SeqIO interface, as
shown below.

The SAM (Sequence Alignment/Map) format and its binary variant, BAM, are a
generic format for storing large nucleotide sequence alignments as generated
in second generation sequencing. They can also hold unaligned reads. The
Bio.SeqIO module in Biopython allows memory efficient access to the reads in
a SAM or BAM file, allowing for example easy conversion to FASTQ:

    >>> from Bio import SeqIO
    >>> SeqIO.convert("SamBam/ex1.sam", "sam", "temp_ex1.fastq", "fastq")
    3270

Or, from the BAM file instead:

    >>> from Bio import SeqIO
    >>> SeqIO.convert("SamBam/ex1.bam", "bam", "temp_ex1.fastq", "fastq")
    3270
    >>> import os
    >>> os.remove("temp_ex1.fastq")

Looking at the first read directly:

    >>> print SeqIO.parse("SamBam/ex1.sam", "sam").next().format("fastq")
    @EAS56_57:6:190:289:82/1
    CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA
    +
    <<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;
    <BLANKLINE>

    >>> print SeqIO.parse("SamBam/ex1.bam", "bam").next().format("fastq")
    @EAS56_57:6:190:289:82/1
    CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA
    +
    <<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;
    <BLANKLINE>

You can also get random access to the reads in a SAM or BAM file by their ID
(but there is no support for indexing by the mapping position):

    >>> from Bio import SeqIO
    >>> sam_reads = SeqIO.index("SamBam/ex1.sam", "sam")
    >>> print sam_reads["EAS56_57:6:190:289:82/2"].format("fastq")
    @EAS56_57:6:190:289:82/2
    AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC
    +
    <<<<<<;<<<<<<<<<<;<<;<<<<;8<6;9;;2;
    <BLANKLINE>

    >>> from Bio import SeqIO
    >>> bam_reads = SeqIO.index("SamBam/ex1.bam", "bam")
    >>> print bam_reads["EAS56_57:6:190:289:82/2"].format("fastq")
    @EAS56_57:6:190:289:82/2
    AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC
    +
    <<<<<<;<<<<<<<<<<;<<;<<<<;8<6;9;;2;
    <BLANKLINE>

For more indepth handling of SAM/BAM files, see the SAM tools project. Their
C API has a Python wrapper called pysam which is very helpful for more
complex task such as accessing columns of the alignment etc:
http://samtools.sourceforge.net
"""

import gzip
import struct

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import SANGER_SCORE_OFFSET

def _decode_flag(flag):
    paired = flag & 0x1
    if paired:
        well_mapped_pair = bool(flag & 0x2)
        mate_unmapped = bool(flag & 0x8)
        mate_strand = bool(flag & 0x20)
        pair1 = bool(flag & 0x40)
        pair2 = bool(flag & 0x80)
    else:
        well_mapped_pair = False
        assert not (flag & 0x2)
        mate_unmapped = False
        assert not flag & 0x8
        mate_rev_strand = False
        assert not flag & 0x20
        pair1 = False
        assert not flag & 0x40
        pair2 = False
        assert not flag & 0x80
        assert not (pair1 and pair2)
    unmapped = bool(flag & 0x4)
    rev_strand = bool(flag & 0x10)
    not_primary = bool(flag & 0x100)
    fails_qc = bool(flag & 0x200)
    duplicate = bool(flag & 0x400)
    return pair1, pair2

def SamIterator(handle, alphabet=generic_dna):
    #Precomputing FASTQ style encoding. Can we share this code with QualityIO?
    q_mapping = dict()
    for letter in range(0, 255):
        q_mapping[chr(letter)] = letter-SANGER_SCORE_OFFSET
    for line in handle:
        if line[0] == "@":
            #Ignore any optional header
            continue
        name, flag, ref_name, ref_pos, mapping_quality, cigar, \
        mate_reference_seq_name, mate_ref_pos, inferred_insert_size, \
        seq_string, quality_string, tags = line.rstrip("\r\n").split("\t",11)
        tags = tags.split("\t")
        
        #If sequence has "." in it, means matches the reference...
        qualities = [q_mapping[letter] for letter in quality_string]

        yield _make_seq_record(name, seq_string, alphabet, qualities, flag)

def _make_seq_record(name, sequence, alphabet, qualities, flag):
    identifier = name
    pair1, pair2 = _decode_flag(int(flag))
    if pair1:
        identifier += "/1"
    elif pair2:
        identifier += "/2"
    return SeqRecord(Seq(sequence, alphabet),
                    id=identifier, name=name, description="",
                    letter_annotations={"phred_quality":qualities})

def _bam_file_header(handle):
    """Read in a BAM file header (PRIVATE).

    Assumings the handle is at the start of the file and has already been
    decompressed (e.g. with the gzip module).

    >>> handle = gzip.open("SamBam/ex1.bam")
    >>> header, num_refs = _bam_file_header(handle)
    >>> header
    ''
    >>> num_refs
    2
    """
    magic = handle.read(4)
    if magic != "BAM\1":
        raise ValueError("After decompression BAM files should start "
                         "with 'BAM\1', not %s" % repr(magic))
    assert 4 == struct.calcsize("<i")
    data = handle.read(4)
    #raise ValueError("Got %s" % repr(data))
    header_length = struct.unpack("<i", data)[0]
    header = handle.read(header_length).rstrip("\0")
    data = handle.read(4)
    #raise ValueError("Got %s" % repr(data))
    num_refs = struct.unpack("<i", data)[0]
    return header, num_refs

def _bam_file_reference(handle):
    """Parse next reference in a BAM file (PRIVATE).

    Assumings the handle is just after the header and has already been                                             
    decompressed (e.g. with the gzip module)

    >>> handle = gzip.open("SamBam/ex1.bam")
    >>> header, num_refs = _bam_file_header(handle)
    >>> for i in range(num_refs):
    ...     print _bam_file_reference(handle)
    ('chr1', 1575)
    ('chr2', 1584)

    """
    ref_name_len = struct.unpack("<i", handle.read(4))[0]
    ref_name = handle.read(ref_name_len).rstrip("\0")
    ref_len = struct.unpack("<i", handle.read(4))[0]
    return ref_name, ref_len

def _build_decoder():
    #TODO - Extend this to cover the ambiguous bases
    answer = {}
    for i,first in enumerate("ACGTN"):
        if first=="N":
            #16+32+64+128 = 240
            value = 240
        else:
            value = 2**(i+4)
        answer[value] = first+"-"
        for j,second in enumerate("ACGTN"):
            if second=="N":
                #1+2+4+8 = 15
                answer[value + 15] = first+second
            else:
                answer[value + 2**j] = first+second
    assert answer[79] == "GN"
    assert answer[240] == "N-"
    assert answer[241] == "NA"
    return answer
_decode_dibase_byte = _build_decoder()


def _bam_file_read_header(handle):
    """Parse the header of the next read in a BAM file (PRIVATE).

    >>> handle = gzip.open("SamBam/ex1.bam")
    >>> header, num_refs = _bam_file_header(handle)
    >>> for i in range(num_refs):
    ...     print _bam_file_reference(handle)
    ('chr1', 1575)
    ('chr2', 1584)
    >>> print _bam_file_read_header(handle)[:3]
    ('EAS56_57:6:190:289:82', 38, 153)
    >>> handle.seek(153)
    >>> print _bam_file_read_header(handle)[:3]
    ('EAS56_57:6:190:289:82', 153, 292)
    >>> handle.seek(153)
    >>> print _bam_file_read_header(handle)[:3]
    ('EAS56_57:6:190:289:82', 153, 292)

    Returns a tuple of the read name, start offset, end offset, etc.

    The offset information is used for indexing - the start offset is to find
    the record on demand, the end offset is used when building the index to
    skip over the rest of the read.
    
    The end offset is also used to determine the end of the tags section.
    """
    #TODO - Check BBH really works for the bin_mq_ml field, defined by
    # bin_mq_nl = bin<<16|mapQual<<8|read_name_len (including NULL)
    start_offset = handle.tell()
    
    fmt = "<iiiBBHHHiiii"
    assert 36 == struct.calcsize(fmt)
    data = handle.read(36)
    if not data:
        raise StopIteration
    if len(data) < 26:
        raise ValueError("Premature end of file")
    #raise ValueError("Data %s = %s" % (repr(data), repr(struct.unpack(fmt, data))))

    block_size, ref_id, ref_pos, read_name_len, bin, map_qual, \
    cigar_len, flag, read_len, mate_ref_id, mate_ref_pos, \
    inferred_insert_size = struct.unpack(fmt, data)

    if read_name_len > 50 or read_name_len <= 0:
        raise ValueError("A read name length of %i probably means the "
                         "parser is out of sync somehow. Read starts:\n%s"
                         "\nand the read name would be %s etc." \
                         % (read_name_len, repr(data), repr(handle.read(25))))
    if read_len > 5000 or read_len <= 0:
        raise ValueError("A read length of %i probably means the parser is out "
                         "of sync somehow. Read starts:\n%s\nand the read name "
                         "would be %s." \
                         % (read_len, repr(data), repr(handle.read(read_name_len))))

    read_name = handle.read(read_name_len).rstrip("\0")
    end_offset = start_offset + block_size + 4
    return read_name, start_offset, end_offset, ref_id, ref_pos, bin, \
           map_qual, cigar_len, flag, read_len, mate_ref_id, mate_ref_pos


def _bam_file_read(handle):
    """Parse the next read in a BAM file (PRIVATE).

    >>> handle = gzip.open("SamBam/ex1.bam")
    >>> header, num_refs = _bam_file_header(handle)
    >>> for i in range(num_refs):
    ...     print _bam_file_reference(handle)
    ('chr1', 1575)
    ('chr2', 1584)
    >>> print _bam_file_read(handle)[:2]
    ('EAS56_57:6:190:289:82', 'CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA')
    >>> print _bam_file_read(handle)[:2]
    ('EAS56_57:6:190:289:82', 'AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC')
    >>> print _bam_file_read(handle)[:2]
    ('EAS51_64:3:190:727:308', 'GGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGG')

    """
    read_name, start_offset, end_offset, ref_id, ref_pos, \
        bin, map_qual, cigar_len, flag, read_len, mate_ref_id, \
        mate_ref_pos = _bam_file_read_header(handle)
    
    fmt = "<%iI" % cigar_len
    cigar = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))
    #raise ValueError("flag %i, cigar length %i, %s" % (flag, cigar_len, repr(cigar)))

    fmt = "<%iB" % ((read_len+1)/2) # round up to make it even
    try:
        binary_seq = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))
    except Exception, err:
        raise ValueError("Expected read name %s of length %i\n%s" \
                         % (repr(read_name), read_len, err))
    seq_string = "".join(_decode_dibase_byte[byte] for byte in binary_seq)
    seq_string = seq_string[:read_len] # remove extra value if odd

    fmt = "<%iB" % read_len
    quals = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))

    #TODO - Parse the tags
    #For now, just skip the tags
    handle.seek(end_offset)

    return read_name, seq_string, quals, flag

def BamIterator(handle, alphabet=generic_dna):
    h = gzip.GzipFile(fileobj=handle)
    header, ref_count = _bam_file_header(h)
    #Skip any reference information
    for i in range(ref_count):
        ref_name, ref_len = _bam_file_reference(h)
    #Loop over the reads
    while True:
        name, seq_string, qualities, flag = _bam_file_read(h)
        yield _make_seq_record(name, seq_string, alphabet, qualities, flag)
    raise StopIteration

def _test():
    """Run the Bio.SeqIO.BamSamIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..", "..", "Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..", "..", "Tests"))
        assert os.path.isfile("SamBam/ex1.sam")
        assert os.path.isfile("SamBam/ex1.bam")
        doctest.testmod(verbose=0)
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
        
if __name__ == "__main__":
    _test()
