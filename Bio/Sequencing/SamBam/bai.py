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
    data = handle.read(4)
    n_ref = struct.unpack("<i", data)[0]
    print "%i references" % n_ref


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
