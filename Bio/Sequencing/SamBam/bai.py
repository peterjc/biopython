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

def _test_bai(handle):
    """Test function for loading a BAI file.

    >>> handle = open("SamBam/ex1.bam.bai", "rb")
    >>> _test_bai(handle)
    2 references
    0 unmapped reads
    >>> handle.close()

    >>> handle = open("SamBam/tags.bam.bai", "rb")
    >>> _test_bai(handle)
    2 references
    0 unmapped reads
    >>> handle.close()

    >>> handle = open("SamBam/bins.bam.bai", "rb")
    >>> _test_bai(handle)
    3 references
    12 unmapped reads
    >>> handle.close()

    """
    indexes, unmapped = _load_bai(handle)
    print "%i references" % len(indexes)
    if unmapped is None:
        print "Index missing unmapped reads count"
    else:
        print "%i unmapped reads" % unmapped

def _load_bai(handle):
    indexes = []
    magic = handle.read(4)
    if magic != "BAI" + chr(1):
        raise ValueError("BAM index files should start 'BAI\1', not %r" \
                         % magic)
    assert 4 == struct.calcsize("<i")
    assert 8 == struct.calcsize("<Q")
    data = handle.read(4)
    n_ref = struct.unpack("<i", data)[0]
    #print "%i references" % n_ref
    for n in xrange(n_ref):
        indexes.append(_load_ref_index(handle))
    #This is missing on very old samtools index files,
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
    """
    #First the chunks for each bin,
    n_bin = struct.unpack("<i", handle.read(4))[0]
    chunks_dict = dict()
    for b in xrange(n_bin):
        bin, chunks = struct.unpack("<ii", handle.read(8))
        chunks_list = []
        for chunk in xrange(chunks):
            #Append tuple of (chunk beginning, chunk end)
            chunks_list.append(struct.unpack("<QQ", handle.read(16)))
        chunks_dict[bin] = chunks_list
    #Now the linear index (for the smallest bins)
    n_intv = struct.unpack("<i", handle.read(4))[0]
    return chunks_dict, struct.unpack("<%iQ" % n_intv, handle.read(8*n_intv))

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
