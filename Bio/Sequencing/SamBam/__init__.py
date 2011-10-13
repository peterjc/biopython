# Copyright 2010-2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
r"""Fairly low level API for working with SAM/BAM files.

This is intended to be written in Pure Python (so that it will work
under PyPy, Jython, etc) but will attempt to follow the pysam API
somewhat (which is a wrapper for the samtools C API).

The low level objects can be used as follows:

>>> data = "frag_5022\t16\tNC_000913_bb\t1\t255\t36M1S\t*\t0\t0\tTCTATTCATTATCTCAATAGCTTTTCATTCTGACTGN\tMMMMMMMMMMMMMMKKKK##%')+.024JMMMMMMM!\tRG:Z:Solexa_test\n"
>>> read = SamRead(data)
>>> read.tid
'frag_5022'
>>> read.flag
16

"""

import gzip
import struct

def SamIterator(handle):
    """Loop over a SAM file returning SamRead objects.

    >>> with open("SamBam/ex1.sam") as handle:
    ...     for read in SamIterator(handle):
    ...         print read.tid, read.flag, read.seq
    ...         if read.tid == "EAS219_FC30151:3:40:1128:1940": break
    EAS56_57:6:190:289:82 69 CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA
    EAS56_57:6:190:289:82 137 AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC
    EAS51_64:3:190:727:308 99 GGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGG
    EAS112_34:7:141:80:875 99 AGCCGAGTCACGGGGTTGCCAGCACAGGGGCTTAA
    EAS219_FC30151:3:40:1128:1940 163 CCGAGTCACGGGGTTGCCAGCACAGGGGCTTAACC

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
    ...         print read.tid, read.flag
    ...         if read.tid == "EAS219_FC30151:3:40:1128:1940": break
    EAS56_57:6:190:289:82 69
    EAS56_57:6:190:289:82 137
    EAS51_64:3:190:727:308 99
    EAS112_34:7:141:80:875 99
    EAS219_FC30151:3:40:1128:1940 163

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
            mate_ref_pos = _bam_file_read_header(h)
        #raw_cigar = h.read(cigar_len * XXX)
        #raw_seq = h.read((read_len+1)/2) # round up to make it even
        #raw_qual = h.read(read_len)
        #TODO - the tags
        h.seek(end_offset)
        yield SamRead("%s\t%i\t...\n" % (read_name, flag))


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
    if len(data) < 26:
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
    return read_name, start_offset, end_offset, ref_id, ref_pos, bin, \
           map_qual, cigar_len, flag, read_len, mate_ref_id, mate_ref_pos



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
        into fields, for later parsing on demand if a property is
        accessed:

        >>> data = 'rd01\t...\n'
        >>> read = SamRead(data)
        >>> print read.tid
        rd01

        Note that a potentially unexpected side effect of this is that
        a malformed entry (e.g. a non-numeric mapping position) may
        not be detected unless accessed:

        >>> print read.flag
        Traceback (most recent call last):
        ...
        ValueError: invalid literal for int() with base 10: '...'

        You can modify values too, this overrides any values parsed
        from the raw data and will be used if saving the record back
        to disk later on:

        >>> read.tid = "Fred"
        >>> print read.tid
        Fred

        """
        self._data = data.rstrip("\n").split("\t")

    def _get_tid(self):
        try:
            return self._tid
        except AttributeError:
            tid = self._data[0]
            self._tid = tid
            return tid
    def _set_tid(self, value):
        self._tid = value
    tid = property(fget = _get_tid, fset = _set_tid,
                   doc = "Template ID (read ID aka QNAME)")

    def _get_flag(self):
        try:
            return self._flag
        except AttributeError:
            flag = int(self._data[1])
            self._flag = flag
            return flag
    def _set_flag(self, value):
        self._flag = value
    flag = property(fget = _get_flag, fset = _set_flag,
                    doc = "FLAG (integer representing bit field)")

    def _get_seq(self):
        try:
            return self._seq
        except AttributeError:
            seq = self._data[9]
            if seq == "*":
                seq = None
            self._seq = seq
            return seq
    def _set_seq(self, value):
        self._seq = value
    seq = property(fget = _get_seq, fset = _set_seq,
                   doc = "SEQ - read sequence bases as string, including soft clipped bases (None if not present)")


#TODO - BamRead class, a subclass where the data is decoded using struct

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
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
    elif os.path.isdir(os.path.join("Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"

if __name__ == "__main__":
    _test()
