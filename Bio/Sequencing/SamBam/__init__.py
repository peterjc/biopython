# Copyright 2010-2012 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
r"""Fairly low level API for working with SAM/BAM files.

This is intended to be written in Pure Python (so that it will work
under PyPy, Jython, etc) but will attempt to follow the pysam API
somewhat (which is a wrapper for the samtools C API).

The SAM and BAM parsers return SamRead and BamRead objects, but these
should behave identically (modulo the tag ordering):

    >>> sam = SamIterator(open("SamBam/ex1.sam"))
    >>> bam = BamIterator(open("SamBam/ex1.bam", "rb"))
    >>> for s, b in zip(sam,bam):
    ...     assert s.qname == b.qname
    ...     assert s.flag == b.flag
    ...     assert s.seq == b.seq
    ...     assert s.qual == b.qual
    ...     assert s.cigar == b.cigar
    ...     assert sorted(s.tags.items()) == sorted(b.tags.items())

Let's look at the FLAG values for the last read, including a couple
of the helper properties to access the bit values (here bit 0x4 is
set for unmapped reads):

    >>> print s.flag, hex(s.flag), s.is_mapped, s.is_unmapped
    83 0x53 True False
    >>> s.is_mapped = False
    >>> print s.flag, hex(s.flag), s.is_mapped, s.is_unmapped
    87 0x57 False True

Here's a sneaky trick, you can access the header like so:

    >>> print repr(sam.text)
    ''

Tricked you - this SAM file has no header! Unsurprisingly, neither does
the BAM version of the file - but the BAM format does separately list
all the reference names and their lengths so a minimal SAM style header
can be inferred automatically:

    >>> print repr(bam.text)
    '@SQ\tSN:chr1\tLN:1575\n@SQ\tSN:chr2\tLN:1584\n'

That is probably atypical though - you would expect a proper SAM header
to be present, indeed it is essential if using things like read groups:

    >>> sam = SamIterator(open("SamBam/tags.sam"))
    >>> print repr(sam.text)
    '@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n'

Or on the BAM equivalent to this SAM file:

    >>> bam = BamIterator(open("SamBam/tags.bam", "rb"))
    >>> print repr(bam.text)
    '@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n'

In addition to the raw header string, you can get the reference names,
and lengths, as tuple, plus the number of references:

    >>> bam.nreferences
    2
    >>> bam.references
    ('chr1', 'chr2')
    >>> bam.lengths
    (100, 200)

The above works on BAM files even if there is no embedded SAM header.
For SAM files it depends on the @SQ lines to work:

    >>> sam.nreferences
    2
    >>> sam.references
    ('chr1', 'chr2')
    >>> sam.lengths
    (100, 200)

The iterators allow you to filter by flag, much like the samtools view
command does. For example, with no filtering we get all the records:

    >>> def count_sam(filename, required_flag=0, excluded_flag=0):
    ...     count = 0
    ...     with open(filename) as handle:
    ...         for read in SamIterator(handle, required_flag, excluded_flag):
    ...             count += 1
    ...     return count
    >>> count_sam("SamBam/ex1.sam")
    3270

The optional argument required_flag selectes only records where those FLAG
bits are set (like how "samtools view -f FLAG ..." with a lower case -f
works). For instance, an argument of 0x2 (hex) or 2 (decimal) requires the
reads be properly mapped:

    >>> count_sam("SamBam/ex1.sam", required_flag=0x2)
    3124

The same applies to the BAM parser as well,

    >>> def count_bam(filename, required_flag=0, excluded_flag=0):
    ...     count = 0
    ...     with open(filename, "rb") as handle:
    ...         for read in BamIterator(handle, required_flag, excluded_flag):
    ...             count += 1
    ...     return count
    >>> count_bam("SamBam/ex1.bam", required_flag=0x2)
    3124

The excluded_flag argument allows you to specify a black list, so to exclude
any properly paired reads (like "samtools view -F FLAG ..." with a capital
-F switch):

   >>> count_sam("SamBam/ex1.sam", excluded_flag=0x2)
   146
   >>> count_bam("SamBam/ex1.bam", excluded_flag=0x2)
   146
   >>> 3270 - 3124
   146

To demonstrate combining flags, consider that flag bit 0x40 indicates
the first read of a pair (or set):

    >>> count_sam("SamBam/ex1.sam", required_flag=0x40)
    1636

If you wanted to select only those reads which are properly mapped (0x2)
and are the first read (0x40), you would use a flag of 0x2 + 0x40 = 0x42,
which in decimal is 2 + 64 = 66.

    >>> count_sam("SamBam/ex1.sam", required_flag=0x42)
    1564

Similarly, the second read in a pair has flag bit 0x80 set, so for just
the read in second properly mapped pairs use 0x2 + 0x80 or 0x82,

    >>> count_sam("SamBam/ex1.sam", required_flag=0x82)
    1560

If both 0x40 and 0x80 are set, the read is neither the first not the last
read in a multi-read fragment (e.g. middle reads in a strobed read).

    >>> count_sam("SamBam/ex1.sam", required_flag=0x120)
    0

You can of course combined the required and excluded flags, for instance
0x4 indicates the read is mapped to the reverse strand, and 0x8 indicates
the read's partner is mapped to the reverse strand. We can use these to
count pairs mapped to the same or different strands. We add the extra
requirement of 0x2 for properly mapped pairs, and 0x40 for first read:

   >>> print "1--> 2-->", count_sam("SamBam/ex1.sam", required_flag=0x42, excluded_flag=0x4+0x8)
   1--> 2--> 1564

   >>> print "1--> <--2", count_sam("SamBam/ex1.sam", required_flag=0x42+0x8, excluded_flag=0x4)
   1--> <--2 0

   >>> print "<--1 2-->", count_sam("SamBam/ex1.sam", required_flag=0x42+0x4, excluded_flag=0x8)
   <--1 2--> 0

   >>> print "<--1 <--2", count_sam("SamBam/ex1.sam", required_flag=0x42+0x4+0x8)
   <--1 <--2 0

Can you guess what kind of paired reads there were? My guess from that
would be 454 paired end reads, since they are on the same strand.

"""

import struct
import sys

#Biopython imports:
from bai import _load_bai
from Bio import bgzf
from Bio._py3k import _as_bytes, _as_string
_empty_bytes_string = _as_bytes("")
_null_byte = _as_bytes("\0")
_ff_byte = _as_bytes("\xFF")

#TODO - Define CIGAR operator constants here?

class SamIterator(object):
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
    >>> print read.rnext, read.pnext #RNEXT and PNEXT aka MRNM and MPOS
    = 290
    >>> print read.tlen #TLEN aka ISIZE
    214
    >>> print read.cigar
    [(0, 35)]
    >>> print read.cigar_str
    35M

    Optional argument required_flag is used like "samtools view -f FLAG ..."
    to only show records where those FLAG bits are set. Thus for example,
    using an argument of 2 requires the reads be properly mapped.

    If your SAM file is gzip compressed, use the Python standard library gzip
    to decompress it during parsing. i.e. handle = gzip.open(filename, "r")

    You should be alerted to mixing up SAM and BAM as inputs,

    >>> sam = SamIterator(open("SamBam/ex1.bam"))
    Traceback (most recent call last):
    ...
    ValueError: Not a SAM file, perhaps it is BAM format or compressed?.

    >>> sam = SamIterator(open("SamBam/ex1.uncompressed.bam"))
    Traceback (most recent call last):
    ...
    ValueError: Not a SAM file, perhaps it is BAM format or compressed?.

    """
    def __init__(self, handle, required_flag=0, excluded_flag=0):
        self._handle = handle
        self._required_flag = required_flag
        self._excluded_flag = excluded_flag
        headers = []
        self._saved_line = None
        self._references = []
        try:
            line = handle.readline()
        except UnicodeDecodeError:
            #This could be almost any binary file, but I want the same error on python 3
            raise ValueError("Not a SAM file, perhaps it is BAM format or compressed?.")
        if sys.version_info[0] < 3 and line.startswith("\x1f\x8b"):
            #This is triggered on Python 2, means BAM or GZIP but want same error as Python 3
            #raise ValueError("This looks like a BAM file or gzipped file, not a SAM file.")
            raise ValueError("Not a SAM file, perhaps it is BAM format or compressed?.")
        if sys.version_info[0] < 3 and line.startswith("BAM" + chr(1)):
            #i.e. a Naked uncompressed BAM file without the gzip/BGZF wrapper
            #raise ValueError("This looks like an uncompressed BAM file, not a SAM file.")
            raise ValueError("Not a SAM file, perhaps it is BAM format or compressed?.")
        while line:
            if line[0] != "@":
                self._saved_line = line
                break
            headers.append(line)
            line = handle.readline()
        self.text = "".join(headers)
        self._references = _sam_header_to_ref_list(self.text)

    @property
    def nreferences(self):
        """Number of reference sequences (read only integer)."""
        return len(self._references)

    @property
    def references(self):
        """Names of the reference sequences (read only tuple)."""
        return tuple(r for r,l in self._references)

    @property
    def lengths(self):
        """Lengths of the reference sequences (read only tuple)."""
        return tuple(l for r,l in self._references)

    @property
    def header(self):
        """Header parsed into two-level dictionary (read only)."""
        return ParseSamHeader(self.text)

    def __iter__(self):
        handle = self._handle
        #Mess about with the first read as a special case since it was
        #already taken from the handle during __iter__ header parsing.
        if self._saved_line:
            line = self._saved_line
            self._saved_line = None
        else:
            line = handle.readline()
        if self._required_flag or self._excluded_flag:
            #Filter the reads
            required_flag = self._required_flag
            excluded_flag = self._excluded_flag
            while line:
                if line[0] == "@":
                    raise ValueError("SAM header @ lines must be before the reads")
                flag = int(line.split("\t",2)[1])
                if (flag & required_flag == required_flag) \
                and not (flag & excluded_flag):
                    yield SamRead(line)
                line = handle.readline()
        else:
            #Take all the reads
            while line:
                if line[0] == "@":
                    raise ValueError("SAM header @ lines must be before the reads")
                yield SamRead(line)
                line = handle.readline()

class BamIterator(object):
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
    >>> print read.rnext, read.pnext #aka RNEXT and PNEXT
    chr1 290
    >>> print read.tlen #aka TLEN
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

    By default a normal BAM file is assumed, using the gzip/BGZF wrapper
    (which may be uncompressed). This parser will also accept 'naked' BAM
    files without the gzip/BGZF wrapping if you use the optional argument
    gzipped=False (such BAM files are very unusual though and is intended
    for debugging low level BAM problems and testing).

    >>> import gzip
    >>> count = 0
    >>> handle = gzip.open("SamBam/ex1.bam", "rb")
    >>> for read in BamIterator(handle, gzipped=False):
    ...     count += 1
    >>> handle.close()
    >>> count
    3270

    """
    def __init__(self, handle, required_flag=0, excluded_flag=0, gzipped=True, bai_filename=None):
        self._handle = handle
        self._required_flag = required_flag
        self._excluded_flag = excluded_flag
        if gzipped:
            h = bgzf.BgzfReader(fileobj=handle, mode="rb")
        else:
            #Uncompressed BAM (useful in testing)
            h = handle
        self.text, ref_count = _bam_file_header(h)
        #Load any reference information
        self._references = [_bam_file_reference(h) for i in range(ref_count)]
        if self.text:
            alt_refs = _sam_header_to_ref_list(self.text)
            if not alt_refs:
                #Append minimal SAM @SQ lines (good idea?)
                self.text += _ref_list_to_sam_header(self._references)
            elif alt_refs != self._references:
                raise ValueError("BAM reference names & lengths inconsistent with SAM header @SQ lines") 
        else:
            #Generate a minimal SAM style header from the BAM header (good idea?)
            self.text = _ref_list_to_sam_header(self._references)
        self._h = h
        if bai_filename:
            #TODO - Make this automatics when we have a normal file handle?
            self._load_index(bai_filename)


    @property
    def nreferences(self):
        """Number of reference sequences (read only integer)."""
        return len(self._references)

    @property
    def references(self):
        """Names of the reference sequences (read only tuple)."""
        return tuple(r for r,l in self._references)

    @property
    def lengths(self):
        """Lengths of the reference sequences (read only tuple)."""
        return tuple(l for r,l in self._references)

    @property
    def header(self):
        """Header parsed into two-level dictionary (read only)."""
        return ParseSamHeader(self.text)

    def __iter__(self):
        h = self._h
        references = self._references
        required_flag = self._required_flag
        excluded_flag = self._excluded_flag
        #Assumes the handle is just after the header!
        #Loop over the reads
        while True:
            read = _load_next_bam_read(h, references, required_flag, excluded_flag)
            #If the read didn't match the required flags will get None
            if read is not None:
                yield read

    def _load_index(self, bai_filename):
        handle = open(bai_filename, "rb")
        indexes, unmapped = _load_bai(handle)
        handle.close()
        if len(indexes) != len(self._references):
            raise ValueError("Have %i references from BAM, but %i from BAI" \
                             % (len(self._references), len(indexes)))
        self._indexes = indexes
        self._unmapped = unmapped

    def fetch(self, reference, start, end):
        """Fetch aligned reads within region specified by reference, start and end.

        Note that start and end should use Python style counting.

        This requires a BAM index file to be present (i.e. a BAI file).
        Currently any BAM index file (BAI file) must be specified explicitly,
        this is an interim measure while indexing is being developed:

        >>> handle = open("SamBam/ex1.bam", "rb")
        >>> bam = BamIterator(handle, bai_filename="SamBam/ex1.bam.bai")
        >>> print bam.nreferences
        2
        >>> print bam.references
        ('chr1', 'chr2')
        >>> print bam._unmapped
        0
        >>> reads = list(bam.fetch("chr1", 120, 150))
        >>> len(reads)
        19
        >>> handle.close()

        """
        #TODO - Nice slice style API using start:end as well/instead?
        references = self._references
        i = self.references.index(reference)
        if not self._indexes:
            raise ValueError("BAI file not loaded")
        index, offsets, mapped, unmapped, u_start, u_end = self._indexes[i]
        #Now find the bins for this region, and the offsets for these bins
        bins = reg2bins(start, end)
        h = self._h
        required_flag = 0
        excluded_flag = 0x4 #Unmapped
        #print "%s region %i:%i covers bins %r" % (reference, start, end, bins)
        #print "Baby-bin offsets:", offsets
        all_chunks = []
        for bin in bins:
            #print "bin %i" % bin
            try:
                chunks = index[bin]
            except KeyError:
                #print "No chunks for bin %i" % bin
                continue
            if bin < 4681:
                min_offset = None
            else:
                #The high number bins 4681-37449 are the baby-bins, 19Kbp,
                #which get a linear index of offsets as well
                min_offset = offsets[bin-4681]
            for s_offset, e_offset in chunks:
                if e_offset < min_offset:
                    assert False
                elif s_offset < min_offset:
                    all_chunks.append((min_offset, e_offset))
                else:
                    all_chunks.append((s_offset, e_offset))
        all_chunks.sort()
        #TODO - Must we remove/merge overlapping chunks?
        for s_offset, e_offset in all_chunks:
                #from Bio.bgzf import split_virtual_offset
                #print "bin %i chunk from virtual offset %i (%i, %i) to %i (%i, %i)" \
                #      % (bin,
                #         s_offset, split_virtual_offset(s_offset)[0], split_virtual_offset(s_offset)[1],
                #         e_offset, split_virtual_offset(e_offset)[0], split_virtual_offset(e_offset)[1])
                if min_offset and s_offset < min_offset:
                    #print "Over-riding with baby-bin's start %i" % min_offset
                    h.seek(min_offset)
                else:
                    h.seek(s_offset)
                #Exploits the fact that virtual offsets are ordered
                while h.tell() < e_offset:
                    #now = h.tell()
                    #print "Now at %i (%i, %i) and about to load read..." \
                    #      % (now, split_virtual_offset(now)[0], split_virtual_offset(now)[1])
                    try:
                        read = _load_next_bam_read(h, references,
                                                   required_flag, excluded_flag)
                    except StopIteration:
                        #Could happen with a bug in BGZF handling.
                        data = h.read(20)
                        now = h.tell()
                        raise ValueError("Premature end of file? Now %i (%i, %i), %r, expected to get to %i (%i, %i)" \
                              % (now, split_virtual_offset(now)[0], split_virtual_offset(now)[1], data,
                                 e_offset, split_virtual_offset(e_offset)[0], split_virtual_offset(e_offset)[1]))
                        break
                    if read is None:
                        #If the read didn't match the required flags will get None
                        pass
                    elif end <= read.pos or read.aend <= start:
                        #Doesn't overlap
                        pass
                    else:
                        assert reference == read.rname, "%s vs %s for read %s" % (reference, read.rname, read.qname)
                        #Looks good
                        yield read
            #print "Done all chunks for bin %i" % bin


def _load_next_bam_read(h, references, required_flag, excluded_flag):
            """Load next BAM read/fragement from handle or raise StopIteration."""
            #TODO - Remove indent (left in for a small diff)
            #TODO - Skip unwanted reads within this function?
            read_name, start_offset, end_offset, ref_id, ref_pos, \
                bin, map_qual, cigar_len, flag, read_len, mate_ref_id, \
                mate_ref_pos, inferred_insert_size, tag_len \
                = _bam_file_read_header(h)
            raw_cigar = h.read(cigar_len * 4)
            raw_seq = h.read((read_len+1)//2) # round up to make it even
            raw_qual = h.read(read_len)
            raw_tags = h.read(tag_len)
            #Can't do offset arithmatic with BGZF
            #assert h.tell() == end_offset, \
            #    "%i vs %i diff %i\n" % (h.tell(), end_offset, h.tell()-end_offset)
            if required_flag and flag & required_flag != required_flag:
                return None
            if flag & excluded_flag:
                return None
            if ref_id == -1:
                ref_name = "*"
            else:
                ref_name = references[ref_id][0]
            if mate_ref_id == -1:
                mate_ref_name = "*"
            else:
                mate_ref_name = references[mate_ref_id][0]
            """
            #Debug code for testing BIN
            if cigar_len:
                #Must decode cigar...
                mapped_len = 0
                for value in struct.unpack("<%iI" % cigar_len, raw_cigar):
                    length = value >> 4
                    value -= length << 4
                    if value in [0,2,3,7,8]:
                        mapped_len += length
                expected_bin = reg2bin(ref_pos, ref_pos+mapped_len)
                if bin != expected_bin:
                    print "%s mapping to %i:%i, expected bin %i, got %i" \
                    % (read_name, ref_pos, ref_pos+mapped_len, expected_bin, bin)
            """
            return BamRead(read_name, flag, ref_name, ref_pos, map_qual,
                          raw_cigar, mate_ref_name, mate_ref_pos,
                          inferred_insert_size, read_len,
                          raw_seq, raw_qual, raw_tags)

def ParseSamHeader(text):
    """Parse the SAM plain text header into a two-level dictionary."""
    d1 = dict()
    for line in text.split("\n"):
        if not line.strip():
            continue
        assert line[0] == "@"
        assert line[3] == "\t"
        k1 = line[1:3]
        d2 = dict()
        for part in line[4:].rstrip().split("\t"):
            assert part[2] == ":"
            k, v = part.split(":",1)
            assert len(k)==2
            if k1=="SQ" and k=="LN":
                #Currently the only case need to cast
                v = int(v)
            d2[k] = v
        #Need to understand when pysam uses a list if want to match...
        if k1 in ["HD"]:
            assert k1 not in d1
            d1[k1] = d2
        else:
            try:
                d1[k1].append(d2)
            except KeyError:
                d1[k1] = [d2]
    return d1

def _ref_list_to_sam_header(references):
    return "".join(["@SQ\tSN:%s\tLN:%i\n" % (name, length) \
                    for name, length in references])

def _sam_header_to_ref_list(header):
    references = []
    for line in header.split("\n"):
        if line.startswith("@SQ\t"):
            r = None
            l = None
            for part in line[3:].rstrip().split("\t"):
                if part.startswith("SN:"):
                    r = part[3:]
                elif part.startswith("LN:"):
                    l = int(part[3:])
            if r is None or l is None:
                raise ValueError("Malformed @SQ header (SN and LN required):\n%r" % line)
            references.append((r,l))
    return references

def reg2bin(beg, end):
    """Turn a beg:end region into a bin BAM/UCSC indexing bin number.

    Based on the C function reg2bin given in the SAM/BAM specification.
    Note that this indexing scheme is limited to references of 512Mbps
    (that is 2^29 base pairs).

    >>> 4681 == reg2bin(9, 13)
    True

    """
    assert 0 <= beg <= end < 2**29, "Bad region %i:%i" % (beg, end)
    end -= 1
    if (beg>>14 == end>>14): return ((1<<15)-1)/7 + (beg>>14)
    if (beg>>17 == end>>17): return ((1<<12)-1)/7 + (beg>>17)
    if (beg>>20 == end>>20): return ((1<<9)-1)/7  + (beg>>20)
    if (beg>>23 == end>>23): return ((1<<6)-1)/7  + (beg>>23)
    if (beg>>26 == end>>26): return ((1<<3)-1)/7  + (beg>>26)
    return 0

def reg2bins(beg, end):
    """Turn beg:end region into list of BAM/UCSC indexing bin numbers overlapping it.

    Based on the C function reg2bins given in the SAM/BAM specification.
    Note that this indexing scheme is limited to references of 512Mbps
    (that is 2^29 base pairs).

    >>> reg2bins(9, 13)
    [0, 1, 9, 73, 585, 4681]

    """
    assert 0 <= beg <= end < 2**29, "Bad region %i:%i" % (beg, end)
    bins = [0]
    end -= 1
    for power, offset in [(26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)]:
        for k in range(offset + (beg>>power), offset + 1 + (end>>power)):
            bins.append(k)
    return bins

    
def _bam_file_header(handle):
    """Read in a BAM file header (PRIVATE).

    Assumings the handle is at the start of the file and has already been
    decompressed (e.g. with the gzip module).

    >>> import gzip
    >>> handle = gzip.open("SamBam/ex1.bam", "rb")
    >>> header, num_refs = _bam_file_header(handle)
    >>> header
    ''
    >>> num_refs
    2
    >>> handle.close()

    """
    magic = handle.read(4)
    if magic != _as_bytes("BAM\1"):
        raise ValueError("After decompression BAM files should start "
                         "with bytes 'BAM\1', not %r" % magic)
    assert 4 == struct.calcsize("<i")
    data = handle.read(4)
    #raise ValueError("Got %s" % repr(data))
    header_length = struct.unpack("<i", data)[0]
    header = _as_string(handle.read(header_length).rstrip(_null_byte))
    data = handle.read(4)
    #raise ValueError("Got %s" % repr(data))
    num_refs = struct.unpack("<i", data)[0]
    return header, num_refs

def _bam_file_reference(handle):
    """Parse next reference in a BAM file (PRIVATE).

    Assumings the handle is just after the header and has already been                                             
    decompressed (e.g. with the gzip module)

    >>> import gzip
    >>> handle = gzip.open("SamBam/ex1.bam", "rb")
    >>> header, num_refs = _bam_file_header(handle)
    >>> for i in range(num_refs):
    ...     print _bam_file_reference(handle)
    ('chr1', 1575)
    ('chr2', 1584)
    >>> handle.close()

    """
    ref_name_len = struct.unpack("<i", handle.read(4))[0]
    ref_name = _as_string(handle.read(ref_name_len).rstrip(_null_byte))
    ref_len = struct.unpack("<i", handle.read(4))[0]
    return ref_name, ref_len

def _bam_file_read_header(handle):
    """Parse the header of the next read in a BAM file (PRIVATE).

    >>> import gzip
    >>> handle = gzip.open("SamBam/ex1.bam", "rb")
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
    >>> handle.close()

    Returns a tuple of the read name, start offset, end offset, etc.

    The offset information is used for indexing - the start offset is to find
    the record on demand, the end offset is used when building the index to
    skip over the rest of the read.
    
    The end offset is also used to determine the end of the tags section.
    """
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

    if read_name_len > 200 or read_name_len <= 0:
        raise ValueError("A read name length of %i probably means the "
                         "parser is out of sync somehow. Read starts:\n%s"
                         "\nand the read name would be %s... with flag %i" \
                         % (read_name_len, repr(data),
                            repr(handle.read(25)), flag))
    #FLAG 516 is for embedding a reference sequence in SAM/BAM
    #FLAG 768 is for a dummy read (usually SEQ is * though)
    if (read_len > 5000 and flag!=516 and flag&768!=768) or read_len < 0:
        raise ValueError("A sequence length of %i probably means the parser is out "
                         "of sync somehow. Read starts:\n%s\nand the read name "
                         "would be %s with flag %i" \
                         % (read_len, repr(data),
                            repr(handle.read(read_name_len)), flag))

    read_name = _as_string(handle.read(read_name_len).rstrip(_null_byte))
    end_offset = start_offset + block_size + 4
    #Block size includes misc fields, read name, seq len, qual len and cigar len
    tag_len = block_size - 32 - read_name_len - ((read_len+1)//2) - read_len - cigar_len * 4
    return read_name, start_offset, end_offset, ref_id, ref_pos, bin, \
           map_qual, cigar_len, flag, read_len, mate_ref_id, mate_ref_pos, \
           inferred_insert_size, tag_len

def _build_decoder():
    decode = {}
    encode = {}
    for i,first in enumerate("=ACMGRSVTWYHKDBN"):
        b = chr(i*16)
        encode[first] = b
        encode[first.lower()] = b
        for j,second in enumerate("=ACMGRSVTWYHKDBN"):
            b = chr(i*16 + j)
            decode[b] = first+second #For Python 2
            decode[i*16+j] = first+second #For Python 3
            encode[first+second] = b
            encode[first.lower()+second] = b
            encode[first.lower()+second.lower()] = b
            encode[first+second.lower()] = b

    assert decode[chr(0)] == "=="
    assert decode[chr(1)] == "=A"
    assert decode[chr(2)] == "=C"
    assert decode[chr(4)] == "=G"
    assert decode[chr(8)] == "=T"
    assert decode[chr(15)] == "=N"

    assert decode[chr(16)] == "A="
    assert decode[chr(32)] == "C="
    assert decode[chr(64)] == "G="
    assert decode[chr(128)] == "T="
    assert decode[chr(240)] == "N="

    assert decode[chr(16+2)] == "AC", decode[chr(16+2)]
    assert decode[chr(128+8)] == "TT", decode[chr(128+8)]
    assert decode[chr(64+4)] == "GG", decode[chr(64+4)]
    assert decode[chr(64+15)] == "GN", decode[chr(64+15)]

    assert encode["gn"] == chr(79)
    assert encode["Gn"] == chr(79)
    assert encode["gN"] == chr(79)
    assert encode["GN"] == chr(79)

    return decode, encode
_decode_dibase_byte, _encode_dibase_byte = _build_decoder()

def _decode_seq(binary, seq_len):
    r"""Helper function to decode BAM style sequence (PRIVATE).

    >>> binary = _as_bytes(chr(16+2) + chr(64+8))
    >>> print _decode_seq(binary, 4)
    ACGT
    >>> print _decode_seq(binary, 3)
    ACG

    TODO - Check last char is equals sign for odd sequences?
    """
    seq = "".join(_decode_dibase_byte[b] for b in binary)
    return seq[:seq_len]
assert _decode_seq(_as_bytes('\x12H'), 4) == 'ACGT'
assert _decode_seq(_as_bytes('\x12H'), 3) == 'ACG'

def _encode_seq(seq):
    r"""Helper function to encode BAM style sequence (PRIVATE).

    >>> _encode_seq("ACGT") == '\x12H'
    True
    >>> _encode_seq("ACG=") == _encode_seq("ACG") == '\x12@'
    True

    """
    answer = []
    for i in range(0, len(seq), 2):
        answer.append(_encode_dibase_byte[seq[i:i+2]])
    return "".join(answer)
assert _encode_seq("ACGT") == '\x12H', _encode_seq("ACGT")
assert _encode_seq("ACG=") == '\x12@', _encode_seq("ACG=")
assert _encode_seq("ACG") == '\x12@'

class SamBamReadTags(dict):
    r"""Represents the tags for a SAM/BAM read.

    Although left ambiguous in early versions of the SAM/BAM specification,
    it is now made clear that each two-letter tag should only appear once.
    Therefore we can use the two-letter tags like dictionary keys.

    >>> tags = SamBamReadTags()
    >>> tags["CO"] = ("Z", "My comment")
    >>> tags["CO"]
    ('Z', 'My comment')
    >>> "CO" in tags
    True

    Trying to access undefined tags raises a KeyError,

    >>> tags["RT"]
    Traceback (most recent call last):
    ...
    KeyError: 'RT'

    Where this first differs from a normal Python dictionary is you cannot
    use anything you like for a key. SAM/BAM tags must be two letters:

    >>> tags["x"] = "Hello"
    Traceback (most recent call last):
    ...
    ValueError: Tag keys must be two letters in SAM/BAM, not 'x'

    Also, the values must always be two-entry tuples, where the first
    component is the SAM/BAM type (e.g. "Z" for strings, "i" for int32,
    and "Bi" for an array of int32 values), and the second component is
    the value (e.g. a string, integer, or list of integers). Basic type
    checking is done:

    >>> tags["xx"] = ("Z", "Hello")
    >>> tags["xx"] = ("A", "H")
    >>> tags["xx"] = ("Bi", [1,2,3])

    Next, if you print the tags or use str(tags), you get a SAM formated
    string, tab separated:

    >>> for tag in sorted(str(tags).split("\t")):
    ...     print tag
    CO:Z:My comment
    xx:B:i,1,2,3

    Note the order of the tags is not important in SAM/BAM.
    """
    def __setitem__(self, key, value):
        if len(key) != 2:
            raise ValueError("Tag keys must be two letters in SAM/BAM, not %r" % key)
        try:
            code, data = value
        except ValueError:
            raise ValueError("Tag value must be 2-tuple of SAM/BAM type and data value")
        if code == "Z":
            if not isinstance(data, basestring):
                raise TypeError("Z type SAM/BAM tag values must be strings, not:\n%r" % data)
        elif code == "A":
            if not isinstance(data, basestring) or len(data)!=1:
                raise TypeError("Z type SAM/BAM tag values must be single character strings, not:\n%r" % data)
        elif code == "i":
            try:
                data = int(data)
            except ValueError:
                raise TypeError("i type SAM/BAM tag values must be integers, not:\n%r" % data)
        elif code == "f":
            try:
                data = float(data)
            except ValueError:
                raise TypeError("f type SAM/BAM tag values must be floats, not:\n%r" % data)
        elif code == "Bf":
            try:
                data = [float(v) for v in data]
            except ValueError:
                raise TypeError("B float arrays require numerical values, not:\n%r" % data)
            except AttributeError:
                raise TypeError("B arrays require a list/tuple/sequence of numbers, not:\n%r" % data)
        elif code.startswith("B"):
            #TODO - bounds checking
            try:
                data = [int(v) for v in data]
            except ValueError:
                raise TypeError("B integer arrays require integer values, not:\n%r" % data)
            except AttributeError:
                raise TypeError("B arrays require a list/tuple/sequence of numbers, not:\n%r" % data)
        elif code=="H":
            try:
                int(data, 16)
            except ValueError:
                raise TypeError("H type SAM/BAM tag values must be hexadecimal strings, not:\n%r" % data)
        else:
            raise ValueError("SAM/BAM tag type %r not supported" % code)
        try:
            order = self._order
        except AttributeError:
            self._order = order = self.keys()
        if key not in order:
            order.append(key)
            assert key in self._order
        return dict.__setitem__(self, key, (code, data))

    def iteritems(self):
        #TODO - Fix this hack, used to preserve tag order for testing
        try:
            order = self._order
        except AttributeError:
            order = []
        for key in order:
            try:
                yield key, self[key]
            except KeyError:
                pass
        for key, value in dict.iteritems(self):
            if key not in order:
                yield key, value

    def __str__(self):
        """Returns the tags tab separated in SAM formatting."""
        tags = []
        for key, (code, data) in self.iteritems():
            if code.startswith("B"):
                tags.append("%s:B:%s,%s" % (key, code[1:], ",".join(str(v) for v in data)))
            else:
                tags.append("%s:%s:%s" % (key, code, data))
        #tags.sort() #Want this consistent regardless of Python implementation
        return "\t".join(tags)

    def _as_bam(self):
        """Returns the tags binary encodes in BAM formatting (bytes string)."""
        data = _empty_bytes_string
        for key, (code, value) in self.iteritems():
            assert len(key) ==2, key
            if code in ["Z", "H"]:
                #Store as null terminated string, easy
                assert isinstance(value, basestring), "%s %s %r" % (key, code, value)
                data += _as_bytes(key + code + value) + _null_byte
            elif code == "i":
                #TODO - Time this with try/except letting struct do bounds checking
                #Integer, but how much space will we need in BAM?
                if value >= 0:
                    #Use an unsigned int
                    if value < 2**8:
                        data += _as_bytes(key + "C") + struct.pack("<B", value)
                    elif value < 2**16:
                        data += _as_bytes(key + "S") + struct.pack("<H", value)
                    elif value < 2**32:
                        data += _as_bytes(key + "I") + struct.pack("<I", value)
                    else:
                        raise ValueError("%s:%s:%i too big for BAM unsigned 32bit int" % (key, code, value))
                else:
                    #Use a signed int
                    if -2**7 < value:
                        data += _as_bytes(key + "c") + struct.pack("<b", value)
                    elif -2**15 < value:
                        data += _as_bytes(key + "s") + struct.pack("<h", value)
                    elif -2**31 < value:
                        data += _as_bytes(key + "i") + struct.pack("<i", value)
                    else:
                        raise ValueError("%s:%s:%i too big for BAM signed 32bit int" % (key, code, value))
            elif code=="f":
                #Float
                data += _as_bytes(key + code) + struct.pack("<f", value)
            elif code[0]=="B":
                #Array of ints/floats etc
                assert len(code)==2
                bam2struct = {"c":"b", "C":"B",
                              "s":"h", "S":"H",
                              "i":"i", "I":"I",
                              "f":"f"}
                data += _as_bytes(key + code) + struct.pack("<I%i%s" % (len(value), bam2struct[code[1]]),
                                                 len(value), *value)
            elif code=="A":
                assert len(value)==1 and isinstance(value, basestring)
                data += _as_bytes(key + code + value)
            else:
                #TODO - BAM encoding of other tag types
                raise NotImplementedError("TODO - Storing %s tags in BAM (here for tag %s)" % (code, key))
        return data


class SamBamRead(object):
    r"""Represents a SAM/BAM entry, i.e. a single read.

    Paired end reads are described by two lines in SAM with the same
    template ID, and become two SamRead objects. Likewise a strobe-read
    or other mulit-read structure becomes multiple SamRead objects.

    This should be API equivalent to the pysam.AlignedRead object.

    This is intended for end users to create read objects from
    scratch, as opposed to the SamRead and BamRead classes which
    are used to create functionally identical classes from raw
    data in SAM and BAM format.

    For example, to represent an unmapped read, we need to use FLAG
    0x4, thus:

    >>> read = SamBamRead(qname="demo", flag=0x4)
    >>> str(read)
    'demo\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n'

    There are helper methods to access the FLAG bit values, e.g.

    >>> print read.flag, hex(read.flag), read.is_unmapped, read.is_qcfail
    4 0x4 True False
    >>> read.is_qcfail = True
    >>> print read.flag, hex(read.flag), read.is_unmapped, read.is_qcfail
    516 0x204 True True
    >>> read.is_qcfail = False
    >>> print read.flag, hex(read.flag), read.is_unmapped, read.is_qcfail
    4 0x4 True False
    >>> read.is_unmapped = False
    >>> print read.flag, hex(read.flag), read.is_unmapped, read.is_qcfail
    0 0x0 False False

    Just printing or using str(read) gives the read formatted as a
    SAM line (with the newline character included).
    """
    def __init__(self, qname="*", flag=0, rname="*", pos=-1,
                 mapq=255, cigar_str="*", rnext="*", pnext=-1,
                 tlen=0, seq=None, qual=None, tags=None):
        #TODO - Store qname, rname, rnext as None rather than *
        #TODO - Type checking
        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.pos = pos
        self.mapq = mapq
        self.cigar_str = cigar_str
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual
        if tags is not None and not isinstance(tags, SamBamReadTags):
            raise TypeError("Bad tags argument")
        self._tags = tags #Only create object on demand

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
                 str(self.mapq), self.cigar_str, self.rnext, str(self.pnext+1),
                 str(self.tlen), seq, qual]
        try:
            #See if this is a SAM record where we never un-parsed the tags
            parts.extend(self._raw_tags)
        except AttributeError:
            if self.tags:
                parts.append(str(self.tags))
        if self.rnext == self.rname and self.rnext != "*":
            parts[6] = "="
        try:
            return "\t".join(parts) + "\n"
        except TypeError, e:
            raise TypeError("%s from join on %r" % (e, parts))

    def _as_bam(self, refs):
        """Returns BAM formatted entry as a bytes string (PRIVATE).

        Required argument 'refs' are the list of reference names in the BAM header.

        TODO? Take lengths argument too and validate mapping positions?
        """
        if self.rname != "*":
            ref_index = refs.index(self.rname)
        else:
            ref_index = -1
        if self.seq:
            l_seq = len(self.seq)
        else:
            assert not self.qual
            #Encoding is MIDNSHP=X
            #Want 'M' = 0, 'I' = 1, 'S' = 4, '=' = 7, 'X' = 8
            #l_seq = sum(op_len for op, op_len in self.cigar if op in [0,1,4,7,8])
            #But for BAM we just want the length for the SEQ and QUAL entries
            l_seq = 0
        if self.pos < 0:
            bin = 0
        elif self.cigar:
            #Encoding is MIDNSHP=X
            #Want 'M' = 0, 'D' = 2, 'N' = 3, '=' = 7, 'X' = 8
            bin = reg2bin(self.pos, self.pos + sum(op_len for op, op_len in self.cigar if op in [0,2,3,7,8]))
        else:
            #Have no CIGAR for unmapped reads (given their partner's position)
            #Do what samtools does:
            bin = reg2bin(self.pos, self.pos + 1)
        bin_mq_nl = int(bin)<<16 | int(self.mapq)<<8 | (len(self.qname)+1)
        try:
            cigar = self._binary_cigar #See BamRead subclass
        except AttributeError:
            if self.cigar:
                cigar = [op_len<<4|op for op, op_len in self.cigar]
                cigar = struct.pack("<%iI" % len(cigar), *cigar)
            else:
                cigar = ""
        flag_nc = self.flag << 16 | (len(cigar)//4)
        if self.rnext == "=":
            next_index = ref_index
        elif self.rnext != "*":
            next_index = refs.index(self.rnext)
        else:
            next_index = -1
        try:
            seq = self._binary_seq #See BamRead subclass
        except AttributeError:
            if self.seq:
                seq = _encode_seq(self.seq)
            else:
                assert l_seq == 0
                assert not self.qual
                seq = ""
        assert len(seq) == (l_seq+1) // 2, "%r (len %i) for %s (len %i = %i)" \
               % (seq, len(seq), self.seq, len(self.seq), l_seq)
        try:
            qual = self._binary_qual #See BamRead subclass
        except AttributeError:
            if self.qual:
                #TODO - Store this as ints? Currently FASTQ encoded...
                #TODO - Reuse the dict in Bio.SeqIO.QualityIO for this
                if self.qual and isinstance(self.qual[0], int):
                    #Iteration over a bytes string gives integers
                    qual = _empty_bytes_string.join(chr(q-33) for q in self.qual)
                else:
                    qual = _as_bytes("".join(chr(ord(q)-33) for q in self.qual))
            else:
                qual = _ff_byte * l_seq
        assert len(qual) == l_seq, "%r (len %i) for %r (len %i = %i)" \
               % (qual, len(qual), self.qual, len(self.qual), l_seq)
        try:
            tags = self._binary_tags #See BamRead subclass
        except AttributeError:
            tags = self.tags._as_bam()
        data = struct.pack("<iiIIiiii",
                           ref_index, self.pos,
                           bin_mq_nl, flag_nc, l_seq,
                           next_index, self.pnext, self.tlen)
        data += _as_bytes(self.qname) + _null_byte + cigar + seq + qual + tags
        #First byte is the length of the REST of the record,
        return struct.pack("<I", len(data)) + data

    #For other FLAG methods, see "magic" after class
    #This one is difference because it flips the bit value,
    #0x4 set means unmapped!
    def _get_mapped(self):
        return not bool(self.flag & 0x4)
    def _set_mapped(self, value):
        if value:
            self.flag = self.flag & ~0x4
        else:
            self.flag = self.flag | 0x4
    is_mapped = property(fget = _get_mapped, fset = _set_mapped,
                         doc = "FLAG 0x200, was this read mapped (bit not set)?")

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

    @property
    def mrnm(self):
        """Use RNEXT instead, replaced mate's reference name (MRNM).

        This is a legacy method for compatibility with early pysam
        (versions 0.5 or older).
        """
        return self.rnext

    @property
    def mpos(self):
        """Use PNEXT instead, replaced mate's position (MPOS).

        This is a legacy method for compatibility with early pysam
        (version 0.5 or older).
        """
        return self.pnext

    @property
    def isize(self):
        """Use template length (TLEN) instead, replaced insert size (ISZIE).

        This is a legacy method for compatibility with early pysam
        (version 0.5 or older).
        """
        return self.tlen

    @property
    def tags(self):
        """Any tags for the read as a dictionary like object."""
        if self._tags:
            return self._tags
        tags = SamBamReadTags()
        self._tags = tags
        return tags

    @property
    def alen(self):
        """Aligned length of the read on the reference genome (None if unavailable).

        The aligned length is given by the sum of the CIGAR M/=/X/D/N operations.

        See also the pos and aend properties.
        """
        if self.cigar is None:
            return None
        length = 0
        for op, op_len in self.cigar:
            #The CIGAR operators MIDNSHP=X are represented using 0 to 7 (as in BAM)
            if op in [0, 2, 3, 6, 7]:
                length += op_len
        return length

    @property
    def aend(self):
        """Aligned end position of the read on the reference (None for unmapped).

        For mapped reads this equals the sum of the pos and alen properties:

        >>> with open("SamBam/ex1.sam") as handle:
        ...     for read in SamIterator(handle):
        ...         if read.is_mapped:
        ...             print read.qname, read.pos, read.alen, read.aend
        ...         if read.qname == "EAS219_FC30151:3:40:1128:1940": break
        EAS56_57:6:190:289:82 99 35 134
        EAS51_64:3:190:727:308 102 35 137
        EAS112_34:7:141:80:875 109 35 144
        EAS219_FC30151:3:40:1128:1940 111 35 146
        """
        start = self.pos
        length = self.alen
        if start is None or length is None:
            return None
        else:
            return start + length

class SamRead(SamBamRead):
    """Represents a SAM/BAM entry created from a SAM file entry.

    This is a subclass of SamBamRead, and it intended for use in creating
    a read object from a line in a SAM file.
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
        >>> print read.pnext
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
        self.rname = parts[2] #RNAME, the actual name, not an integer!
        self.pos = int(parts[3]) - 1
        self.mapq = int(parts[4])
        self.cigar_str = parts[5]
        self.rnext = parts[6] #RNEXT aka MRNM, not an integer!
        self.pnext = int(parts[7]) - 1 #PNEXT aka MPOS
        self.tlen = int(parts[8]) #TLEN aka ISAIZE
        if parts[9] == "*":
            self.seq = None
        else:
            self.seq = parts[9]
        if parts[10] == "*":
            self.qual = None
        else:
            self.qual = parts[10]
        self._raw_tags = parts[11:]
        self._tags = None

    @property
    def tags(self):
        if self._tags:
            return self._tags
        tags = SamBamReadTags()
        for tag in self._raw_tags:
            key, code, data = tag.split(":",2)
            if code=="B":
                assert data[1]==","
                code = code + data[0]
                data = data[2:].split(",")
            #The tags object should do any casting for us
            tags[key] = (code, data)
        self._tags = tags
        del self._raw_tags
        return tags


class BamRead(SamBamRead):
    """Represents a SAM/BAM entry created from a BAM file entry.

    This is a subclass of SamBamRead, and it intended for use in creating
    a read object from an entry in a BAM file.
    """
    def __init__(self, qname, flag, rname, pos, mapq, binary_cigar,
                 rnext, pnext, tlen, read_len, binary_seq, binary_qual,
                 binary_tags):
        r"""Create a BamRead object.

        This is a lazy-parsing approach to loading SAM/BAM files, so
        all the parser does is grab the raw data and pass it to this
        object. The bare minimum parsing is done at this point.

        >>> read = BamRead(qname='rd01', flag=1, rname="*", pos=-1, mapq=255, binary_cigar='', rnext="*", pnext=-1, tlen=0, read_len=0, binary_seq='', binary_qual='', binary_tags='')
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
        self.rname = rname #RNAME, the actual name, not the integer!
        self.pos = pos
        self.mapq = mapq
        self._binary_cigar = binary_cigar
        self.rnext = rnext #RNEXT aka MRNM, the actual name, not the integer!
        self.pnext = pnext #PNEXT aka MPOS
        self.tlen = tlen #TLEN aka ISIZE
        self._read_len = read_len
        self._binary_seq = binary_seq
        self._binary_qual = binary_qual
        self._binary_tags = binary_tags
        self._tags = None
    
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
            seq = _decode_seq(self._binary_seq, self._read_len)
            self._seq = seq
            return seq
    def _set_seq(self, value):
        del self._binary_seq
        self._seq = value
    seq = property(fget = _get_seq, fset = _set_seq,
                   doc = "SEQ - read sequence bases as string, including soft clipped bases")

    def _get_qual(self):
        try:
            return self._qual
        except AttributeError:
            #TODO - Reuse dict mapping from FASTQ parser, will be faster
            if sys.version_info[0] >= 3:
                #Iteration over a bytes string gives integers
                qual = "".join(chr(33+byte) for byte in self._binary_qual)
            else:
                qual = "".join(chr(33+ord(byte)) for byte in self._binary_qual)
            self._qual = qual
            return qual
    def _set_qual(self, value):
        del self._binary_qual
        self._qual = value
    qual = property(fget = _get_qual, fset = _set_qual,
                    doc = "QUAL - read quality as FASTQ encoded string, including soft clipped bases")

    @property
    def tags(self):
        if self._tags:
            return self._tags
        raw = self._binary_tags
        tags = SamBamReadTags()
        while raw:
            tag, code, value, raw = _next_tag_raw(raw)
            tags[tag] = (code, value)
        self._tags = tags
        del self._binary_tags
        return tags

#Magic to define lots of very similar properties at once:
#(can we do this to the base class in a way that covers the subclasses?)
def _make_prop(bit, prop, help):
    #import sys
    #sys.stderr.write("Making %s property\n" % prop)
    def _get_flag_bit(self):
        return bool(self.flag & bit)
    def _set_flag_bit(self,value):
            if value:
                self.flag = self.flag |bit
            else:
                self.flag = self.flag &~bit
    return property(fget=_get_flag_bit, fset=_set_flag_bit, doc=help)
for bit, prop, help in [
    (0x4, "unmapped", "FLAG 0x4, is this read unmapped?"),
    (0x10, "reverse", "FLAG 0x10, is this read mapped to the reverse strand?"),
    (0x200, "qcfail", "FLAG 0x200, did the read fail quality control (QC)?"),
    (0x400, "duplicate", "FLAG 0x400, is this read a PCR or optical duplicate?"),
    ]:
    p = _make_prop(bit, prop, help)
    setattr(SamBamRead, "is_%s" % prop, p)
    setattr(SamRead, "is_%s" % prop, p)
    setattr(BamRead, "is_%s" % prop, p)
del bit, prop, help, p, _make_prop


def _next_tag_raw(raw):
    tag = _as_string(raw[0:2])
    code = _as_string(raw[2:3])
    if code == "B":
        sub_code = _as_string(raw[3:4])
        length = struct.unpack("<I", raw[4:8])[0]
        if sub_code == "f": #float
            values = struct.unpack(("<%if" % length), raw[8:8+length*4])
            return tag, "Bf", values, raw[8+length*4:]
        elif sub_code == "i": #int32
            values = struct.unpack(("<%ii" % length), raw[8:8+length*4])
            return tag, "Bi", values, raw[8+length*4:]
        elif sub_code == "I": #int32
            values = struct.unpack(("<%iI" % length), raw[8:8+length*4])
            return tag, "BI", values, raw[8+length*4:]
        elif sub_code == "s": #int16
            values = struct.unpack(("<%ih" % length), raw[8:8+length*2])
            return tag, "Bs", values, raw[8+length*2:]
        elif sub_code == "S": #int16
            values = struct.unpack(("<%iH" % length), raw[8:8+length*2])
            return tag, "BS", values, raw[8+length*2:]
        elif sub_code == "c": #int8
            values = struct.unpack(("<%ib" % length), raw[8:8+length])
            return tag, "Bc", values, raw[8+length:]
        elif sub_code == "C": #int8
            values = struct.unpack(("<%iB" % length), raw[8:8+length])
            return tag, "BC", values, raw[8+length:]
        else:
            raise ValueError("Unknown BAM tag B sub-element type %r (for %r tag)" % (sub_code, tag))
    elif code == "C": #u_int8
        return tag, "i", ord(raw[3:4]), raw[4:]
    elif code == "S": #u_int16
        return tag, "i", struct.unpack("<H", raw[3:5])[0], raw[5:]
    elif code == "I": #u_int32
        return tag, "i", struct.unpack("<I", raw[3:7])[0], raw[7:]
    elif code == "s": #int16
        return tag, "i", struct.unpack("<h", raw[3:5])[0], raw[5:]
    elif code == "i": #int32
        return tag, "i", struct.unpack("<i", raw[3:7])[0], raw[7:]
    elif code == "c": #int8
        value = ord(raw[3])
        if value >= 128:
            #Negative bit set
            value -= 256
        return tag, "i", value, raw[4:]
    elif code == "A": #Single char
        return tag, "A", _as_string(raw[3:4]), raw[4:]
    elif code == "Z": #Null terminated string
        i = 3
        while ord(raw[i:i+1]) != 0: i+= 1
        return tag, "Z", _as_string(raw[3:i]), raw[i+1:]
    elif code == "H": #Hex, null terminated string
        i = 3
        while ord(raw[i:i+1]) != 0: i+= 1
        if (i-3) % 2 == 1:
            #Warning only?
            raise ValueError("Odd number of bytes for hex string? %r" % raw)
        return tag, "H", _as_string(raw[3:i]), raw[i+1:]
    elif code == "f": #Single precision float
        #TODO, leave it as a float rather than turning it into a string
        #which is a short term solution during testing
        return tag, "f", str(struct.unpack("<f",raw[3:7])[0]), raw[7:]
    else:
        raise ValueError("Unknown BAM tag element type %r (for %r tag)" % (code, tag))


def SamWriter(handle, reads, header="", referencenames=None, referencelengths=None):
    """Writes a complete SAM file including any header, returns read count.

    >>> bam = BamIterator(open("SamBam/ex1.bam", "rb"))
    >>> handle = open("saved.sam", "w")
    >>> count = SamWriter(handle, bam)
    >>> handle.close()
    >>> print "Saved %i reads to SAM file" % count
    Saved 3270 reads to SAM file

    For a more complicated example, to get only the properly mapped paired
    reads (i.e. where FLAG 0x2 is set):

    >>> bam = BamIterator(open("SamBam/ex1.bam", "rb"), required_flag=0x2)
    >>> handle = open("saved.sam", "w")
    >>> count = SamWriter(handle, bam)
    >>> handle.close()
    >>> print "Saved %i reads to SAM file" % count
    Saved 3124 reads to SAM file

    """
    #Consistency check:
    header, references = _cross_check_header_refs(reads, header,
                                                  referencenames,
                                                  referencelengths)
    #Write SAM header, perform some basic sanity testing:
    handle.write(_check_header_text(header))
    #Write SAM reads:
    count = 0
    for read in reads:
        handle.write(str(read))
        count += 1
    return count

def BamWriter(handle, reads, header="", referencenames=None, referencelengths=None, gzipped=True):
    """Writes a complete BAM file including any header, returns read count.

    Note that if you do not supply any header information, it will be copied
    from the reads object if that is a SamIterator or BamIterator. If it is
    something else like a list or a generator that won't work.

    The header information is optional, but only in the special case of no
    mapped reads. Even then you may want a header for things like read groups.

    >>> sam = SamIterator(open("SamBam/ex1_header.sam"))
    >>> handle = open("saved_from_sam.bam", "wb")
    >>> count = BamWriter(handle, sam)
    >>> handle.close()
    >>> print "Saved %i reads to BAM file" % count
    Saved 3270 reads to BAM file

    And a (silly) example for testing, BAM to BAM,

    >>> bam = BamIterator(open("SamBam/ex1_header.bam", "rb"))
    >>> handle = open("saved_from_bam.bam", "wb")
    >>> count = BamWriter(handle, bam)
    >>> handle.close()
    >>> print "Saved %i reads to BAM file" % count
    Saved 3270 reads to BAM file

    """
    if gzipped:
        handle = bgzf.BgzfWriter(fileobj=handle)
    #Consistency check:
    header, references = _cross_check_header_refs(reads, header,
                                                  referencenames,
                                                  referencelengths)
    header = _check_header_text(header)
    #Write BAM header:
    handle.write("BAM" + chr(1))
    handle.write(struct.pack("<I", len(header)))
    handle.write(header)
    handle.write(struct.pack("<I", len(references)))
    for r, l in references:
        try:
            handle.write(struct.pack("<I", len(r)+1))
            handle.write(r + chr(0))
            handle.write(struct.pack("<I", l))
        except Exception:
            raise ValueError("Problem with reference %r, %r" % (r,l))
    #Want to give it its own BGZF block, even if there is no SAM header:
    handle.flush()
    #Write reads:
    count = 0
    refs = [r for (r,l) in references]
    for read in reads:
        handle.write(read._as_bam(refs))
        count += 1
    handle.flush() #Important with BGZF as it buffers data
    handle.flush() #Again to trigger empty BGZF block as EOF marker
    return count

def _check_header_text(text):
    """Removes blank lines, ensures trailing newlines, spots some errors."""
    lines = []
    for line in text.split("\n"):
        if not line:
            continue
            #raise ValueError("Blank line in header")
        elif line[0] != "@":
            raise ValueError("SAM header lines must start with @, not %s" % line)
        lines.append(line + "\n")
    return "".join(lines)

def _cross_check_header_refs(reads, header="", referencenames=None, referencelengths=None):
    if not header:
        #If the reads argument is a SamIterator or BamIterator this works:
        try:
            header = reads.text
        except AttributeError:
            pass
    if referencenames is not None or referencelengths is not None:
        assert len(referencenames) == len(referencelengths)
        references = zip(referencenames, referencelengths)
        alt_refs = _sam_header_to_ref_list(header)
        if not alt_refs:
            #Append minimal @SQ lines to the header
            header += _ref_list_to_sam_header(references)
        elif alt_refs!=references:
            raise ValueError("Reference names and lengths inconsistent with header @SQ lines")
    else:
        references = _sam_header_to_ref_list(header)
    return header, references


def _comp_float(a,b):
    return repr(a)==repr(b) or a==b or abs(float(a)-float(b))<0.000001,

def _test_misc():
    print "Misc tests..."
    import gzip
    import os
    import sys

    def compare(a_iter, b_iter):
        from itertools import izip_longest
        assert a_iter.nreferences == b_iter.nreferences, \
            "%r vs %r" % (a_iter.nreferences, b_iter.nreferences)
        assert a_iter.references == b_iter.references, \
            "%r vs %r" % (a_iter.references, b_iter.references)
        assert a_iter.lengths == b_iter.lengths, \
            "%r vs %r" % (a_iter.lengths, b_iter.lengths)
        if a_iter.header and b_iter.header:
            assert a_iter.header == b_iter.header, \
                "%r vs %r" % (a_iter.header, b_iter.header)
        for a, b in izip_longest(a_iter, b_iter):
            assert b is not None, "Extra read in a: %r" % str(a)
            assert a is not None, "Extra read in b: %r" % str(b)
            assert a.qname == b.qname, "%r vs %r" % (a.qname, b.qname)
            assert a.flag == b.flag, "%r vs %r" % (a.flag, b.flag)
            assert a.rname == b.rname, "%r vs %r" % (a.rname, b.rname)
            assert a.pos == b.pos, "%r vs %r" % (a.pos, b.pos)
            assert a.mapq == b.mapq, "%r vs %r" % (a.mapq, b.mapq)
            assert a.cigar == b.cigar, "%r vs %r" % (a.cigar, b.cigar)
            assert a.mrnm == b.mrnm, "%r vs %r" % (a.mrnm, b.mrnm)                                                        
            assert a.mpos == b.mpos, "%r vs %r" % (a.pos, b.pos)
            assert a.isize == b.isize, "%r vs %r" % (a.isize, b.isize)
            assert a.seq == b.seq, "%r vs %r" % (a.seq, b.seq)
            assert a.qual == b.qual, "%r vs %r" % (a.qual, b.qual)
            assert a.tags.keys() == b.tags.keys(), "%r vs %r" % (a.tags.keys(), b.tags.keys())
            for key in a.tags:
                assert a.tags[key][0] == b.tags[key][0], \
                       "%s tag %s has codes %s vs %s" % \
                       (a.qname, key, a.tags[key][0], b.tags[key][0])
                if a.tags[key][0]=="f":
                    assert str(a.tags[key][1]) == str(b.tags[key][1]) or \
                           abs(a.tags[key][1] - b.tags[key][1]) < 0.00001, \
                           "%s tag %s have %f vs %f" % \
                           (a.qname, key, a.tags[key][1], b.tags[key][1])
                elif a.tags[key][0]=="Bf":
                    pass
                else:
                    assert a.tags[key][1] == b.tags[key][1], \
                           "%s tag %s:%s has %s vs %s" % \
                           (a.qname, key, a.tags[key][0], a.tags[key][1], b.tags[key][1])
            #Float formating in tags is annoying...
            #assert str(a) == str(b), "Reads disagree,\n%s\n%s\n" % (a,b)
    compare(SamIterator(open("SamBam/tags.sam")),
            BamIterator(open("SamBam/tags.bam", "rb")))
    compare(SamIterator(open("SamBam/tags.sam")),
            BamIterator(gzip.open("SamBam/tags.bam"), gzipped=False))
    for read in SamIterator(open("SamBam/tags.sam")):
        #TODO - Test API for getting tag values
        tag = str(read).rstrip("\n").split("\t")[-1]
        assert read.qname == "tag_" + tag, \
               "%s vs tag of %s" % (read.qname, tag)
    for read in BamIterator(open("SamBam/tags.bam", "rb")):
        tag = str(read).rstrip("\n").split("\t")[-1]
        assert read.qname.startswith("tag_" + tag[:5])
        if read.qname == "tag_" + tag:
            continue
        if ":f:" in tag:
            old = float(read.qname.split(":")[2])
            new = float(tag.split(":")[2])
            assert _comp_float(old, new), \
                   "%s vs tag of %s" % (read.qname, tag)
        elif ":B:" in tag:
            assert read.qname.startswith("tag_" + tag[:5]), \
                   "%s vs tag of %s" % (read.qname, tag)
            if ":B:f," in read.qname:
                old = tuple(float(v) for v in read.qname.split(":B:f,")[1].split(","))
                new = tuple(float(v) for v in tag.split(":")[2][2:].split(","))
                assert len(old) == len(new), \
                       "Count mismatch %i vs %i in %s\n%s\n%s\n" \
                       % (read.qname, len(old), len(new), old, new)
                for a,b in zip(old, new):
                    assert _comp_float(a, b), \
                           "Mismatch in %s,\n%s\n%s" % (read.qname, old, new)
            else:
                #Integers
                old = [int(v) for v in read.qname.split(":B:")[1][2:].split(",")]
                new = [int(v) for v in tag.split(":")[2][2:].split(",")]
                assert old == new, "Mismatch in read %s vs %s\n%r\n%r" % (read.qname, tag, old, new)
        else:
            assert read.qname == "tag_" + tag, \
                   "Please check %s vs tag of %s" % (read.qname, tag)

    try:
        #This is in Python 2.6+, but we need it on Python 3
        from io import BytesIO
    except ImportError:
        #Must be on Python 2.5 or older
        from StringIO import StringIO as BytesIO
    for sam_filename, bam_filename in [
        ("SamBam/tags.sam", "SamBam/tags.bam"),
        ("SamBam/ex1_header.sam", "SamBam/ex1_header.bam"),
        ("SamBam/bins.sam", "SamBam/bins.bam"),
        ]:
        #BAM -> BAM
        if not os.path.isfile(bam_filename):
            continue
        print "%s -> %s check..." % (bam_filename, bam_filename)
        handle = BytesIO()
        count = BamWriter(handle, BamIterator(open(bam_filename, "rb")),
                          gzipped=False)
        if handle.getvalue() != gzip.open(bam_filename).read():
            sys.stderr.write("ERROR: Couldn't reproduce %s -> %s\n" \
                             % (bam_filename, bam_filename))
        #SAM -> BAM
        if not os.path.isfile(sam_filename):
            continue
        print "%s -> %s check..." % (sam_filename, bam_filename)
        handle = BytesIO()
        count = BamWriter(handle, SamIterator(open(sam_filename)),
                          gzipped=False)
        if handle.getvalue() != gzip.open(bam_filename).read():
            sys.stderr.write("ERROR: Couldn't reproduce %s -> %s\n" \
                             % (sam_filename, bam_filename))
        #Comparing SAM output is tricky due to float formatting...

    print "Done"

def _test():
    """Run the module's doctests (PRIVATE).

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..", "..", "..", "Tests")):
        #print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..", "..", "..", "Tests"))
        doctest.testmod()
        print "Done"
        _test_misc()
        os.chdir(cur_dir)
        del cur_dir
    elif os.path.isdir(os.path.join("Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        print "Done"
        _test_misc()
        os.chdir(cur_dir)
        del cur_dir

if __name__ == "__main__":
    _test()
