# Copyright 2010-2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
r"""Fairly low level API for working with SAM/BAM files.

This is intended to be written in Pure Python (so that it will work
under PyPy, Jython, etc) but will attempt to follow the pysam API
somewhat (which is a wrapper for the samtools C API).

>>> from itertools import izip
>>> sam = SamIterator(open("SamBam/ex1.sam"))
>>> bam = BamIterator(open("SamBam/ex1.bam", "rb"))
>>> for s, b in izip(sam,bam):
...     assert s.qname == b.qname
...     assert s.seq == b.seq
...     assert s.qual == b.qual
...     assert s.cigar == b.cigar
...     b.mrnm = s.mrnm #hack for testing
...     assert str(s) == str(b), "SAM:\n%r\nBAM:\n%r\n" % (str(s), str(b))

"""

import gzip
import struct

def SamIterator(handle):
    """Loop over a SAM file returning SamRead objects.

    >>> with open("SamBam/ex1.sam") as handle:
    ...     for read in SamIterator(handle):
    ...         print read.qname, read.flag, read.pos, read.seq
    ...         if read.qname == "EAS219_FC30151:3:40:1128:1940": break
    EAS56_57:6:190:289:82 69 99 CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA
    EAS56_57:6:190:289:82 137 99 AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC
    EAS51_64:3:190:727:308 99 102 GGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGG
    EAS112_34:7:141:80:875 99 109 AGCCGAGTCACGGGGTTGCCAGCACAGGGGCTTAA
    EAS219_FC30151:3:40:1128:1940 163 111 CCGAGTCACGGGGTTGCCAGCACAGGGGCTTAACC
    >>> print read.qual
    <<<<<<<<<<<<<<<<<<<;<<5;;<<<9;;;;7:
    >>> print read.mapq
    99
    >>> print read.rname, read.pos
    chr1 111
    >>> print read.mrnm, read.mpos #aka RNEXT and PNEXT
    = 290
    >>> print read.isize #aka TLEN
    214
    >>> print read.cigar
    [(0, 35)]
    >>> print read.cigar_str
    35M


    >>> count = 0
    >>> with open("SamBam/ex1.sam") as handle:
    ...     for read in SamIterator(handle):
    ...         count += 1
    >>> count
    3270

    """
    for line in handle:
        if line[0] == "@":
            #Ignore any optional header
            continue
        else:
            yield SamRead(line)


def BamIterator(handle):
    """Loop over a BAM file returning BamRead objects.

    This should be functionally identical to using a SAM file
    with the SamIterator function.

    >>> with open("SamBam/ex1.bam", "rb") as handle:
    ...     for read in BamIterator(handle):
    ...         print read.qname, read.flag, read.pos, read.seq
    ...         if read.qname == "EAS219_FC30151:3:40:1128:1940": break
    EAS56_57:6:190:289:82 69 99 CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA
    EAS56_57:6:190:289:82 137 99 AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC
    EAS51_64:3:190:727:308 99 102 GGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGG
    EAS112_34:7:141:80:875 99 109 AGCCGAGTCACGGGGTTGCCAGCACAGGGGCTTAA
    EAS219_FC30151:3:40:1128:1940 163 111 CCGAGTCACGGGGTTGCCAGCACAGGGGCTTAACC
    >>> print read.qual
    <<<<<<<<<<<<<<<<<<<;<<5;;<<<9;;;;7:
    >>> print read.mapq
    99
    >>> print read.rname, read.pos
    chr1 111
    >>> print read.mrnm, read.mpos #aka RNEXT and PNEXT
    chr1 290
    >>> print read.isize #aka TLEN
    214
    >>> print read.cigar
    [(0, 35)]
    >>> print read.cigar_str
    35M


    >>> count = 0
    >>> with open("SamBam/ex1.bam", "rb") as handle:
    ...     for read in BamIterator(handle):
    ...         count += 1
    >>> count
    3270

    """
    h = gzip.GzipFile(fileobj=handle)
    header, ref_count = _bam_file_header(h)
    #Load any reference information
    references = [_bam_file_reference(h) for i in range(ref_count)]
    #Loop over the reads
    while True:
        read_name, start_offset, end_offset, ref_id, ref_pos, \
            bin, map_qual, cigar_len, flag, read_len, mate_ref_id, \
            mate_ref_pos, inferred_insert_size, tag_len \
            = _bam_file_read_header(h)
        raw_cigar = h.read(cigar_len * 4)
        raw_seq = h.read((read_len+1)/2) # round up to make it even
        raw_qual = h.read(read_len)
        raw_tags = h.read(tag_len)
        assert h.tell() == end_offset, \
            "%i vs %i diff %i\n" % (h.tell(), end_offset, h.tell()-end_offset)
        ref_name = references[ref_id][0]
        mate_ref_name = references[mate_ref_id][0]
        yield BamRead(read_name, flag, ref_name, ref_pos, map_qual,
                      raw_cigar, mate_ref_name, mate_ref_pos,
                      inferred_insert_size, read_len,
                      raw_seq, raw_qual, raw_tags)


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
    >>> x = handle.seek(153)
    >>> print _bam_file_read_header(handle)[:3]
    ('EAS56_57:6:190:289:82', 153, 292)
    >>> x = handle.seek(153)
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
    if len(data) < 36:
        raise ValueError("Premature end of file")
    #raise ValueError("Data %s = %s" % (repr(data), repr(struct.unpack(fmt, data))))

    block_size, ref_id, ref_pos, read_name_len, map_qual, bin, \
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
    #Block size includes misc fields, read name, seq len, qual len and cigar len
    tag_len = block_size - 32 - read_name_len - ((read_len+1)/2) - read_len - cigar_len * 4
    return read_name, start_offset, end_offset, ref_id, ref_pos, bin, \
           map_qual, cigar_len, flag, read_len, mate_ref_id, mate_ref_pos, \
           inferred_insert_size, tag_len

def _build_decoder():
    answer = {}
    for i,first in enumerate("=ACMGRSVTWYHKDBN"):
        answer[chr(i*16)] = first+"="
        for j,second in enumerate("=ACMGRSVTWYHKDBN"):
            answer[chr(i*16 + j)] = first+second
    assert answer[chr(79)] == "GN"
    assert answer[chr(240)] == "N="
    assert answer[chr(241)] == "NA"
    return answer
_decode_dibase_byte = _build_decoder()


#TODO - Have a Flag class?

class SamRead(object):
    """Represents a SAM/BAM entry, i.e. a single read.

    Paired end reads are described by two lines in SAM with the same
    template ID, and become two SamRead objects. Likewise a strobe-read
    or other mulit-read structure becomes multiple SamRead objects.

    This should be API equivalent to the pysam.AlignedRead object.
    """
    def __init__(self, data):
        r"""Create a SamRead object.

        This is a lazy-parsing approach to loading SAM/BAM files, so
        all the parser does is grab the raw data and pass it to this
        object. The bare minimum parsing is done - splitting the text
        into fields, tags for instance are parsed on demand if accessed.

        >>> data = "frag_5022\t16\tNC_000913_bb\t1\t255\t36M1S\t*\t0\t0\tTCTATTCATTATCTCAATAGCTTTTCATTCTGACTGN\tMMMMMMMMMMMMMMKKKK##%')+.024JMMMMMMM!\tRG:Z:Solexa_test\n"
        >>> read = SamRead(data)
        >>> print read.qname
        frag_5022
        >>> print read.pos
        0
        >>> print read.mpos
        -1

        Getting back to the SAM formatted line is easy, print it or use str(),

        >>> str(read)
        "frag_5022\t16\tNC_000913_bb\t1\t255\t36M1S\t*\t0\t0\tTCTATTCATTATCTCAATAGCTTTTCATTCTGACTGN\tMMMMMMMMMMMMMMKKKK##%')+.024JMMMMMMM!\tRG:Z:Solexa_test\n"
        >>> data == str(read)
        True

        Note that a potentially unexpected side effect of this is that
        a malformed entry (e.g. invalid tags) may not be detected unless
        accessed.

        You can modify values too, this overrides any values parsed
        from the raw data and will be used if saving the record back
        to disk later on:

        >>> read.qname = "Fred"
        >>> print read.qname
        Fred

        """
        parts = data.rstrip("\n").split("\t")
        self.qname = parts[0]
        self.flag = int(parts[1])
        self.rname = parts[2] #not an integer!
        self.pos = int(parts[3]) - 1
        self.mapq = int(parts[4])
        self.cigar_str = parts[5]
        self.mrnm = parts[6]
        self.mpos = int(parts[7]) - 1 #aka PNEXT
        self.isize = int(parts[8]) #aka TLEN
        if parts[9] == "*":
            self.seq = None
        else:
            self.seq = parts[9]
        if parts[10] == "*":
            self.qual = None
        else:
            self.qual = parts[10]
        self._tags = parts[11:]
    
    def __repr__(self):
        """Returns simple representation of the read for debugging."""
        return "%s(qname=%r, flag=%i, ...)" \
                         % (self.__class__.__name__, self.qname, self.flag)

    def __str__(self):
        """Returns the read as a SAM line (including trailing new line)."""
        seq = self.seq
        if seq is None:
            seq = "*"
        qual = self.qual
        if qual is None:
            qual = "*"
        parts = [self.qname, str(self.flag), self.rname, str(self.pos+1),
                 str(self.mapq), self.cigar_str, self.mrnm, str(self.mpos+1),
                 str(self.isize), seq, qual] + self._tags
        if self.mrnm == self.rname and self.mrnm != "*":
            parts[6] = "="
        try:
            return "\t".join(parts) + "\n"
        except TypeError, e:
            raise TypeError("%s from join on %r" % (e, parts))

    @property
    def cigar(self):
        """CIGAR string parsed into a list of tuples (operator code, count).
        
        The CIGAR operators MIDNSHP=X are represented using 0 to 7 (as in BAM),
        so a cigar string of 36M2I3M becomes [(0, 36), (1, 2), (0, 3)] etc.
        
        Any empty CIGAR string (represented as * in SAM) is given as None.
        """
        cigar = self.cigar_str
        if cigar == "*":
            return None
        answer = []
        count = ""
        for letter in cigar:
            if letter.isdigit():
                count += letter #string addition
            else:
                operator = "MIDNSHP=X".find(letter)
                if operator == -1:
                    raise ValueError("Invalid character %s in CIGAR %s" \
                                     % (letter, cigar))
                answer.append((operator, int(count)))
                count = ""
        return answer

 
class BamRead(SamRead):
    def __init__(self, qname, flag, rname, pos, mapq, binary_cigar, mrnm, mpos, isize, read_len, binary_seq, binary_qual, binary_tags):
        r"""Create a BamRead object.

        This is a lazy-parsing approach to loading SAM/BAM files, so
        all the parser does is grab the raw data and pass it to this
        object. The bare minimum parsing is done at this point.

        >>> read = BamRead(qname='rd01', flag=1, rname="*", pos=-1, mapq=255, binary_cigar='', mrnm="*", mpos=-1, isize=0, read_len=0, binary_seq='', binary_qual='', binary_tags='')
        >>> print read.qname
        rd01

        Note that a potentially unexpected side effect of this is that
        a malformed entry (e.g. a non-numeric mapping position) may
        not be detected unless accessed:

        TODO - example

        Getting the read as a SAM formatted line is easy, print it or use str(),

        >>> str(read)
        'rd01\t1\t*\t0\t255\t*\t*\t0\t0\t\t\n'
    
        You can modify values too, this overrides any values parsed
        from the raw data and will be used if saving the record back
        to disk later on:

        >>> read.qname = "Fred"
        >>> print read.qname
        Fred

        """
        self.qname = qname
        self.flag = flag
        self.rname = rname #the actual name, not the integer!
        self.pos = pos
        self.mapq = mapq
        self._binary_cigar = binary_cigar
        self.mrnm = mrnm #aka RNEXT, the actual name, not the integer!
        self.mpos = mpos #aka PNEXT
        self.isize = isize #aka TLEN
        self._read_len = read_len
        self._binary_seq = binary_seq
        self._binary_qual = binary_qual
        self._binary_tags = binary_tags
    
    @property
    def cigar(self):
        """CIGAR string parsed into a list of tuples (operator code, count).
        
        The CIGAR operators MIDNSHP=X are represented using 0 to 7 (as in BAM).

        Any empty CIGAR string (represented as * in SAM) is given as None.
        """
        cigar = self._binary_cigar
        length = len(cigar) // 4
        if not length:
            return None
        answer = []
        for value in struct.unpack("<%iI" % length, cigar):
            length = value >> 4
            value -= length << 4
            answer.append((value, length))
        return answer

    @property
    def cigar_str(self):
        """CIGAR string as used in the SAM format (string)."""
        cigar = self._binary_cigar
        length = len(cigar) // 4
        if not length:
            return "*"
        answer = ""
        for value in struct.unpack("<%iI" % length, cigar):
            length = value >> 4
            value -= length << 4
            answer += "%i%s" % (length, "MIDNSHP=X"[value])
        return answer

    def _get_seq(self):
        try:
            return self._seq
        except AttributeError:
            seq = "".join(_decode_dibase_byte[byte] for byte in self._binary_seq)
            seq = seq[:self._read_len] # remove extra value if odd
            self._seq = seq
            return seq
    def _set_seq(self, value):
        self._seq = value
    seq = property(fget = _get_seq, fset = _set_seq,
                   doc = "SEQ - read sequence bases as string, including soft clipped bases")

    def _get_qual(self):
        try:
            return self._qual
        except AttributeError:
            #TODO - Reuse dict mapping from FASTQ parser, will be faster
            qual = "".join(chr(33+ord(byte)) for byte in self._binary_qual)
            self._qual = qual
            return qual
    def _set_qual(self, value):
        self._qual = value
    qual = property(fget = _get_qual, fset = _set_qual,
                    doc = "QUAL - read quality as FASTQ encoded string, including soft clipped bases")

    @property
    def _tags(self):
        """Decord the binary tags into a SAM style string (PRIVATE."""
        raw = self._binary_tags
        answer = []
        while raw:
            tag, bam_code, sam_code, value, raw = _next_tag_raw(raw)
            answer.append("%s:%s:%s" % (tag, sam_code, value))
        return answer

def _next_tag_raw(raw):
    tag = raw[0:2]
    code = raw[2]
    if code == "B":
        sub_code = raw[3]
        length = struct.unpack("<I", raw[3:7])
        if sub_code == "S":
            value = raw[7:7+length*unit]
            return tag, value, "Z", raw[7+length*unit:]
        else:
            raise ValueError("Unknown BAM tag B sub-element type %r (for %r tag)" % (sub_code, tag))
    elif code == "C": #u_int8
        return tag, code, "i", ord(raw[3]), raw[4:]
    else:
        raise ValueError("Unknown BAM tag element type %r (for %r tag)" % (code, tag))



def _pysam():
    try:
        import pysam
    except ImportError:
        pass
    print "Running tests against pysam..."
    
    def compare(a_iter, b_iter):
        from itertools import izip
        for a, b in izip(a_iter, b_iter):
             assert a.qname == b.qname, "%r vs %r" % (a.qname, b.qname)
             assert a.flag == b.flag, "%r vs %r" % (a.flag, b.flag)
             assert a.mapq == b.mapq, "%r vs %r" % (a.mapq, b.mapq)
             assert a.cigar == b.cigar, "%r vs %r" % (a.cigar, b.cigar)
             assert a.seq == b.seq, "%r vs %r" % (a.seq, b.seq)
             assert a.qual == b.qual, "%r vs %r" % (a.qual, b.qual)
             #Would compare other fields and str(a)==str(b) but pysam broken,
             #See http://code.google.com/p/pysam/issues/detail?id=74
             #and http://code.google.com/p/pysam/issues/detail?id=75
    
    #TODO, use pysam on the SAM file
    #http://code.google.com/p/pysam/issues/detail?id=73
    #compare(SamIterator(open("SamBam/ex1.sam")),
    #        pysam.Samfile("SamBam/ex1.sam", "r"))
    compare(SamIterator(open("SamBam/ex1.sam")),
            pysam.Samfile("SamBam/ex1.bam", "rb"))
    compare(BamIterator(open("SamBam/ex1.bam", "rb")),
            pysam.Samfile("SamBam/ex1.bam", "rb"))
    print "Done"

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
        _pysam()
        os.chdir(cur_dir)
        del cur_dir
    elif os.path.isdir(os.path.join("Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        print "Done"
        _pysam()
        os.chdir(cur_dir)
        del cur_dir

if __name__ == "__main__":
    _test()
