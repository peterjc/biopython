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
    >>> handle.close()

    >>> handle = open("SamBam/tags.bam.bai", "rb")
    >>> _test_bai(handle)
    2 references
    >>> handle.close()

    >>> handle = open("SamBam/bins.bam.bai", "rb")
    >>> _test_bai(handle)
    3 references
    >>> handle.close()
    """
    magic = handle.read(4)
    if magic != "BAI" + chr(1):
        raise ValueError("BAM index files should start 'BAI\1', not %r" \
                         % magic)
    assert 4 == struct.calcsize("<i")
    assert 8 == struct.calcsize("<Q")
    data = handle.read(4)
    n_ref = struct.unpack("<i", data)[0]
    print "%i references" % n_ref
    for n in xrange(n_ref):
        data = handle.read(4)
        n_bin = struct.unpack("<i", data)[0]
        print " - ref %i has %i bins" % (n, n_bin)
        for b in xrange(n_bin):
            data = handle.read(8)
            bin, chunks = struct.unpack("<ii", data)
            #print "   - bin %i aka %i has %i chunks" % (b, bin, chunks)
            for chunk in xrange(chunks):
                data = handle.read(16)
                chunk_beg, chunk_end = struct.unpack("<QQ", data)
                #print "     - chunk %i from %i to %i" \
                #      %  (chunk, chunk_beg, chunk_end)
        data = handle.read(4)
        n_intv = struct.unpack("<i", data)[0]
        print "    - bin %i aka %i has %i 16kbp intervals" \
              % (b, bin, n_intv)
        data = handle.read(8*n_intv)
        ioffsets = struct.unpack("<%iQ" % n_intv, data)
        #print "      %r" % (ioffsets,)
    #This is missing on very old samtools index files,
    data = handle.read(8)
    if data:
        unmapped = struct.unpack("<Q", data)[0]
        print "%i unmapped reads" % unmapped
    else:
        print "Index missing unmapped reads count"
    data = handle.read()
    if data:
        print "%i extra bytes" % len(data)
        print repr(data)


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
