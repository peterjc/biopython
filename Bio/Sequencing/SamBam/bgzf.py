# Copyright 2010-2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
r"""Fairly low level API for working with BGZF files (e.g. BAM files).

The SAM/BAM file format (Sequence Alignment/Map) comes in a plain text
format (SAM), and a compressed binary format (BAM). The latter uses a
modified form of gzip compression called BGZF, which in principle can
be applied to any file format. This is described together with the
SAM/BAM file format at http://samtools.sourceforge.net/SAM1.pdf


Aim of this module
------------------

The Python gzip library can be used to read BGZF files, since for
decompression they are just (specialised) gzip files. What this
module aims to facilitate is random access to BGZF files (using the
'virtual offset' idea), and writing BGZF files (which means using
suitably sized gzip blocks and writing the extra 'BC' field in the
gzip headers). The existing gzip library will be used internally.

Initially this will be used to provide random access and writing of
BAM files. However, the BGZF format could also be used on other
sequential data (in the sense of one record after another), such
as most of the sequence data formats supported in Bio.SeqIO (like
FASTA, FASTQ, GenBank, etc).


Technical Introduction to BGZF
------------------------------

The gzip file format allows multiple compressed blocks, each of which
could be a stand alone gzip file. As an interesting bonus, this means
you can use Unix "cat" to combined to gzip files into one by
concatenating them. Also, each block can have one of several compression
levels (including uncompressed, which actually takes up a little bit
more space due to the gzip header).

What the BAM designers realised was that random access to data stored
in traditional gzip files was slow, breaking the file into gzip blocks
would allow fast random access to each block. To access a particular
piece of the decompressed data, you just need to know which block it
starts in (the offset of the gzip block start), and how far into the
(decompressed) contents of the block you need to read.

One problem with this is finding the gzip block sizes efficiently.
You can do it with a standard gzip file, but it requires every block
to be decompressed -- and that would be rather slow.

All that differs in BGZF is that compressed size of each gzip block
is limited to 2^16 bytes, and an extra 'BC' field in the gzip header
records this size. Traditional decompression tools can ignore this,
and unzip the file just like any other gzip file.

The point of this is you can look at the first BGZF block, find out
how big it is from this 'BC' header, and thus seek immediately to
the second block, and so on.

The BAM indexing scheme records read positions using a 64 bit
'virtual offset', comprising coffset<<16|uoffset, where uoffset is
the file offset of the BGZF block containing the start of the read
(unsigned integer using up to 64-16 = 48 bits), and coffset is the
offset within the (decompressed) block (unsigned 16 bit integer).

This limits you to BAM files where the last block starts by 2^48
bytes, or 256 petabytes, and the decompressed size of each block
is at most 2^16 bytes, or 64kb. Note that this matches the BGZF
'BC' field size which limits the compressed size of each block to
2^16 bytes, allowing for BAM files to use BGZF with no gzip
compression (useful for intermediate files in memory to reduced
CPU load).

"""
#TODO - Move somewhere else in Bio.* namespace?

import gzip

def open(filename, mode):
    if "r" in mode.lower():
        raise NotImplementedError("TODO - Use gzip.open() for now")
    elif "w" in mode.lower() or "a" in mode.lower():
        return BgzfWriter(filename, mode)
    else:
        raise ValueError("Bad mode %r" % mode)

def make_virtual_offset(block_start_offset, within_block_offset):
    """Compute a BGZF virtual offset from block start and within block offsets.

    The BAM indexing scheme records read positions using a 64 bit
    'virtual offset', comprising in C terms:

    within_block_offset<<16 | block_start_offset

    Here block_start_offset is the file offset of the BGZF block
    start (unsigned integer using up to 64-16 = 48 bits), and
    within_block_offset within the (decompressed) block (unsigned
    16 bit integer).

    >>> make_virtual_offset(0,0)
    0
    >>> make_virtual_offset(0,1)
    1
    >>> make_virtual_offset(0, 2**16 - 1)
    65535
    >>> make_virtual_offset(0, 2**16)
    Traceback (most recent call last):
    ...
    ValueError: Require 0 <= within_block_offset < 2**16

    >>> make_virtual_offset(1,0)
    65536
    >>> make_virtual_offset(1,1)
    65537
    >>> make_virtual_offset(1, 2**16 - 1)
    131071

    >>> make_virtual_offset(100000,0)
    6553600000
    >>> make_virtual_offset(100000,1)
    6553600001
    >>> make_virtual_offset(100000,10)
    6553600010

    >>> make_virtual_offset(2**48,0)
    Traceback (most recent call last):
    ...
    ValueError: Require 0 <= block_start_offset < 2**48

    """
    if within_block_offset < 0 or within_block_offset >= 2**16:
        raise ValueError("Require 0 <= within_block_offset < 2**16")
    if block_start_offset < 0 or block_start_offset >= 2**28:
        raise ValueError("Require 0 <= block_start_offset < 2**48")
    return (block_start_offset<<16) | within_block_offset

def split_virtual_offset(virtual_offset):
    """Divides a 64-bit BGZF virtual offset into block start & within block offsets.

    >>> split_virtual_offset(6553600000)
    (100000, 0)
    >>> split_virtual_offset(6553600010)
    (100000, 10)

    """
    start = virtual_offset>>16
    return start, virtual_offset ^ (start<<16)

class _BgzfWriter(object):
    def __init__(filename, mode):
        if "w" not in mode.lower() \
        and "a" not in mode.lower():
            raise ValueError("Must use write or append mode, not %r" % mode)
        handle = open(filename, mode)
        self._handle = handle
        self._buffer = "" #Bytes string

    def write(self, data):
        pass


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
        del cur_dir
    elif os.path.isdir(os.path.join("Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        print "Done"
        del cur_dir

if __name__ == "__main__":
    _test()
