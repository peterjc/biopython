# Copyright 2012 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
r"""Fairly low level API for working with SAM/BAM index files (BAI).

This is intended to be written in Pure Python (so that it will work
under PyPy, Jython, etc) but will attempt to follow the pysam API
somewhat (which is a wrapper for the samtools C API).
"""

import struct
import gzip

from Bio._py3k import _as_bytes

_BAM_MAX_BIN =  37450 # (8^6-1)/7+1
_BAI_magic = _as_bytes("BAI\1")


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


def _test_bai(handle):
    """Test function for loading a BAI file.

    >>> handle = open("SamBam/ex1.bam.bai", "rb")
    >>> _test_bai(handle)
    2 references
    1 bins, 1 linear baby-bins, 1446 reads mapped, 18 unmapped
    1 bins, 1 linear baby-bins, 1789 reads mapped, 17 unmapped
    0 unmapped unplaced reads
    >>> handle.close()

    >>> handle = open("SamBam/tags.bam.bai", "rb")
    >>> _test_bai(handle)
    2 references
    1 bins, 1 linear baby-bins, 417 reads mapped, 0 unmapped
    0 bins, 0 linear baby-bins, ? reads mapped, ? unmapped
    0 unmapped unplaced reads
    >>> handle.close()

    >>> handle = open("SamBam/bins.bam.bai", "rb")
    >>> _test_bai(handle)
    3 references
    1 bins, 1 linear baby-bins, 194 reads mapped, 42 unmapped
    5 bins, 4 linear baby-bins, 684 reads mapped, 287 unmapped
    73 bins, 64 linear baby-bins, 10196 reads mapped, 5043 unmapped
    12 unmapped unplaced reads
    >>> handle.close()

    """
    indexes, unmapped_unplaced = _load_bai(handle)
    print "%i references" % len(indexes)
    for chunks, linear, mapped, ref_unmapped, u_start, u_end, e_start, e_end, in indexes:
        if mapped is None:
            assert ref_unmapped is None
            print "%i bins, %i linear baby-bins, ? reads mapped, ? unmapped" \
                  % (len(chunks), len(linear))
        else:
            print "%i bins, %i linear baby-bins, %i reads mapped, %i unmapped" \
                  % (len(chunks), len(linear), mapped, ref_unmapped)
        if e_start or e_end:
            print "Embedded reference offset %i to %i" % (e_start, e_end)
    if unmapped_unplaced is None:
        print "Index missing unmapped unplaced reads count"
    else:
        print "%i unmapped unplaced reads" % unmapped_unplaced

def _load_bai(handle):
    indexes = []
    magic = handle.read(4)
    if magic != _BAI_magic:
        raise ValueError("BAM index files should start %r, not %r" \
                         % (_BAI_magic, magic))
    assert 4 == struct.calcsize("<i")
    assert 8 == struct.calcsize("<Q")
    data = handle.read(4)
    n_ref = struct.unpack("<i", data)[0]
    #print "%i references" % n_ref
    for n in xrange(n_ref):
        indexes.append(_load_ref_index(handle))
    #This is missing on very old samtools index files,
    #and isn't in the SAM/BAM specifiction yet either.
    #This was reverse engineered vs "samtools idxstats"
    data = handle.read(8)
    if data:
        unmapped = struct.unpack("<Q", data)[0]
        #print "%i unmapped reads" % unmapped
    else:
        unmapped = None
        #print "Index missing unmapped reads count"
    data = handle.read()
    if data:
        print "%i extra bytes" % len(data)
        print repr(data)
    return indexes, unmapped

def _load_ref_index(handle):
    """Load offset chunks for bins (dict), and linear index (tuple).

    This assumes the handle is positioned at the start of the next
    reference block in the BAI file.

    It returns a dictionary for the chunks, and a list for the linear
    index. The chunk dictionary keys are bin numbers, and its values
    are lists of chunk begining and end virtual offsets. The linear
    index is just a tuple of virtual offsets (the position of the first
    aligned read in that interval) for the smallest sized bins.

    In addition, you also get six values back: The number of mapped
    reads, unmapped reads, the start and end offsets of the reads,
    and the start and end offset for the embedded reference. None
    is used if the information is missing (i.e. out of date index
    file lacking these values).
    """
    mapped = None
    unmapped = None
    unmapped_start = None
    unmapped_end = None
    embed_start = None
    embed_end = None
    #First the chunks for each bin,
    n_bin = struct.unpack("<i", handle.read(4))[0]
    chunks_dict = dict()
    for b in xrange(n_bin):
        bin, chunks = struct.unpack("<ii", handle.read(8))
        if bin == _BAM_MAX_BIN:
            #At the time of writing this isn't in the SAM/BAM specification,
            #gleaned from the samtools source code instead.
            assert chunks >= 2, chunks
            unmapped_start, unmapped_end = struct.unpack("<QQ", handle.read(16))
            mapped, unmapped = struct.unpack("<QQ", handle.read(16))
            if chunks >= 3:
                embed_start, embed_end = struct.unpack("<QQ", handle.read(16))
                assert chunks == 3, chunks
        else:
            chunks_list = []
            for chunk in xrange(chunks):
                #Append tuple of (chunk beginning, chunk end)
                chunks_list.append(struct.unpack("<QQ", handle.read(16)))
            chunks_dict[bin] = chunks_list
    #Now the linear index (for the smallest bins)
    n_intv = struct.unpack("<i", handle.read(4))[0]
    return chunks_dict, struct.unpack("<%iQ" % n_intv, handle.read(8*n_intv)), \
           mapped, unmapped, unmapped_start, unmapped_end, embed_start, embed_end


def idxstats(bam_filename, bai_filename):
    """Generator function returning tuples for each reference.

    Mimics the output of 'samtools idxstats example.bam' returning
    reference name (string), reference length, number of mapped
    reads, and number of placed but unmapped reads (integers).
    Finally returns a tuple of "*", 0, 0, and the number of
    unplaced unmapped reads.

    >>> for values in idxstats("SamBam/ex1.bam", "SamBam/ex1.bam.bai"):
    ...     print "%s\t%i\t%i\t%i" % values
    chr1   1575      1446      18
    chr2   1584      1789      17
    *   0      0      0

    """
    #Don't need random access, so can just use gzip not bgzf
    handle = gzip.open(bam_filename, "rb")
    from Bio.Sequencing.SamBam import _bam_file_header #lazy import
    from Bio.Sequencing.SamBam import _bam_file_reference #lazy import
    header_text, ref_count = _bam_file_header(handle)
    references = [_bam_file_reference(handle) for i in range(ref_count)]
    handle.close()

    handle = open(bai_filename, "rb")
    indexes, unmapped = _load_bai(handle)
    if unmapped is None:
        raise ValueError("Old index lacks unmapped read information, re-index your BAM file")

    if len(indexes) != len(references):
        raise ValueError("BAM file has %i references, BAI has %i" \
                  % (len(references), len(indexes)))

    for (reference, length), (chunks, linear, mapped, ref_unmapped, u_start, u_end, e_start, e_end) in zip(references, indexes):
        if mapped is None:
            mapped = 0
        if ref_unmapped is None:
            ref_unmapped = 0
        if linear:
            #Get one linear index chunk per 16kb (2**14 bp)
            min_len = (2**14) * (len(linear)-1)
            max_len = (2**14) * len(linear)
            if not (min_len <= length <= max_len):
                import warnings
                warnings.warn("WARNING: BAM file says %s is %i bp, but BAI says %i to %i bp"
                              " (from %i linear index entries each of 16384bp)\n" \
                              % (reference, length, min_len, max_len, len(linear)))
        yield reference, length, mapped, ref_unmapped
    yield "*", 0, 0, unmapped


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
