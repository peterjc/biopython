# Copyright 2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
r"""Fairly low level API for working with bzip2 files

This builds on the bz2 library included with Python but also exposes
the block structure to allow more efficient random access. Note that
under Python 2 that using bz2.BZ2File(...) would only decompress the
first block, fixed for Python 3, see http://bugs.python.org/issue1625

The aim of this library is to read bzip2 files and provide a file like
object where seek and tell are modified to exploit the block structure
(see also the Bio.bgzf module).

Warning about namespaces
------------------------

It is considered a bad idea to use "from XXX import *" in Python, because
it pollutes the names space. This is a real issue with Bio.bzip2 (and the
standard Python library gzip) because they contain a function called open
i.e. Suppose you do this:

>>> from Bio.bzip2 import *
>>> print open.__module__
Bio.bzip2

Or,

>>> from gzip import *
>>> print open.__module__
gzip

Notice that the open function has been replaced. You can "fix" this if you
need to by importing the built-in open function:

>>> from __builtin__ import open
>>> print open.__module__
__builtin__

However, what we recommend instead is to keep the name space, e.g.

>>> from Bio import bzip2
>>> print bzip2.open.__module__
Bio.bzip2
>>> print open.__module__
__builtin__

Example
-------

This is an ordinary FASTQ file compressed using a sngle bzip2 block,
so it be decompressed using bz2,

>>> import bz2
>>> handle = bz2.BZ2File("Quality/example.fastq.bz2", "rb")
>>> handle.tell()
0
>>> data = handle.read()
>>> len(data)
234
>>> handle.tell()
234
>>> handle.close()

The next file actually uses multiple bzip2 blocks, which means that
on Python 2 using bz2.BZ2File you'll only get the first block back.
The reader in this module will work on both Python 2 and Python 3.

So far nothing special. But rather than mimicking the traditional
offsets we can ask for and use block-based offsets as tuples.
"""
#TODO - Move somewhere else in Bio.* namespace?

import bz2
import struct
import __builtin__ #to access the usual open function

from Bio._py3k import _as_bytes, _as_string

#For Python 2 can just use: _bgzf_magic = '\x1f\x8b\x08\x04'
#but need to use bytes on Python 3
_bzip2_magic = _as_bytes("BZh")
_empty_bytes_string = _as_bytes("")
_bytes_newline = _as_bytes("\n")

def open(filename, mode="rb"):
    """Open a bzip2 file for reading, writing or appending."""
    if "r" in mode.lower():
        return BZip2Reader(filename, mode)
    elif "w" in mode.lower() or "a" in mode.lower():
        return BZip2Writer(filename, mode)
    else:
        raise ValueError("Bad mode %r" % mode)

def BZip2Blocks(handle):
    """Low level debugging function to inspect bzip2 blocks.

    Returns the block start offset, the block length (add these for
    the start of the next block), and the decompressed length of the
    blocks contents.

    >>> from __builtin__ import open
    >>> handle = open("Quality/example.fastq.bz2", "rb")
    >>> for values in BZip2Blocks(handle):
    ...     print "Raw start %i, raw length %i; data start %i, data length %i" % values
    Raw start 0, raw length 154; data start 0, data length 234
    >>> handle.close()

    The above small file uses one block, and manages to compress
    the 234 bytes of data into just 154 bytes. The next example is
    the same data deliberately split into three blocks - which
    imposes quite a file size overhead:

    >>> handle = open("Quality/example.fastq.3b.bz2", "rb")
    >>> for values in BZip2Blocks(handle):
    ...     print "Raw start %i, raw length %i; data start %i, data length %i" % values
    Raw start 0, raw length 92; data start 0, data length 78
    Raw start 92, raw length 92; data start 78, data length 78
    Raw start 184, raw length 95; data start 156, data length 78
    >>> handle.close()

    """
    chunk = 10000
    raw_data = _empty_bytes_string
    data_start = 0
    while True:
        raw_data += handle.read(chunk - len(raw_data))
        if not raw_data:
            raise StopIteration
        start_offset = handle.tell() - len(raw_data)
        d = bz2.BZ2Decompressor()
        block_length = len(raw_data)
        data_len = len(d.decompress(raw_data))
        while True:
            raw_data = handle.read(chunk)
            if not raw_data:
                break
            block_length += len(raw_data)
            try:
                data += d.decompress(raw_data)
            except bz2.EOFError:
                break
        raw_data = d.unused_data
        block_length -= len(raw_data)
        yield start_offset, block_length, data_start, data_len
        data_start += data_len


def _load_bzip2_block(handle, text_mode=False):
    #Change indentation later...
    magic = handle.read(3)
    if not magic:
        #End of file
        raise StopIteration
    if magic != _bzip2_magic:
        raise ValueError(r"A bzip2 block should start with "
                         r"%r, not %r; handle.tell() now says %r"
                         % (_bzip2_magic, magic, handle.tell()))
    raise NotImplementedError


class BZip2Reader(object):
    r"""BZIP2 reader, acts like a read only handle but seek/tell differ.
    """

    def __init__(self, filename=None, mode=None, fileobj=None, max_cache=3):
        if fileobj:
            assert filename is None and mode is None
            handle = fileobj
            assert "b" in handle.mode.lower()
            self._text = False
        else:
            if "w" in mode.lower() \
            or "a" in mode.lower():
                raise ValueError("Must use read mode (default), not write or append mode")
            self._text = "b" not in mode.lower()
            handle = __builtin__.open(filename, "rb")
        self._handle = handle
        self.max_cache = max_cache
        self._block_start_offset = None
        self._load_block()

    def _load_block(self, start_offset=None):
        raise NotImplementedError

    def readline(self):
        if self._text:
            newline = "\n"
        else:
            newline = _bytes_newline
        i = self._buffer.find(newline, self._within_block_offset)
        if i != -1:
            data = self._buffer[self._within_block_offset:i+1]
            self._within_block_offset = i + 1
            assert data.endswith(newline)
            return data
        else:
            data = self._buffer[self._within_block_offset:]
            self._load_block() #will reset offsets
            if not self._buffer:
                return data #EOF
            else:
                #TODO - Avoid recursion
                return data + self.readline()

    def close(self):
        self._handle.close()
        self._buffer = None
        self._block_start_offset = None


class BgzfWriter(object):

    def __init__(self, filename=None, mode=None, fileobj=None, blocksize=9):
        if fileobj:
            assert filename is None and mode is None            
            handle = fileobj
        else:
            if "w" not in mode.lower() \
            and "a" not in mode.lower():
                raise ValueError("Must use write or append mode, not %r" % mode)
            handle = __builtin__.open(filename, mode)
        self._handle = handle
        self._buffer = _empty_bytes_string
        self.blocksize = blocksize

    def _write_block(self, block):
        raise NotImplementedError

    def close(self):
        if self._buffer:
            self.flush()
        self._handle.close()


def _test():
    import os
    import doctest
    if os.path.isdir(os.path.join("..","Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","Tests"))
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
