# Copyright 2010-2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Fairly low level API for working with SAM/BAM files.

This is intended to be written in Pure Python (so that it will work
under PyPy, Jython, etc) but will attempt to follow the pysam API
somewhat (which is a wrapper for the samtools C API).

"""

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
        not be detected unless accessed.

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
