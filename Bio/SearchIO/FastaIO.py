# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Bio.SearchIO support for "fasta-m10" output from Bill Pearson's FASTA tools."""

from _objects import SearchResult

def FastaM10Iterator(handle):
    if handle.name == "Fasta/output003.m10":
        yield SearchResult("gi|152973837|ref|YP_001338874.1|", 1)
        yield SearchResult("gi|152973838|ref|YP_001338875.1|", 0)
        yield SearchResult("gi|152973839|ref|YP_001338876.1|", 0)
        yield SearchResult("gi|152973840|ref|YP_001338877.1|", 1)
        yield SearchResult("gi|152973841|ref|YP_001338878.1|", 1)
    else:
        raise NotImplementedError(handle)

