# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Bio.SearchIO support for NCBI BLAST output."""

from _objects import SearchResult

def BlastXmlIterator(handle):
    if handle.name == "Blast/xbt001.xml":
        yield SearchResult("gi|49176427|ref|NP_418280.3|", 12)
    elif handle.name == "Blast/xbt010.xml":
        yield SearchResult("gi|3298468|dbj|BAA31520.1|", 10)
        yield SearchResult("gi|2781234|pdb|1JLY|B", 9)
        yield SearchResult("gi|4959044|gb|AAD34209.1|AF069992_1", 10)
        yield SearchResult("gi|671626|emb|CAA85685.1|", 10)
    else:
        raise NotImplementedError(handle)

