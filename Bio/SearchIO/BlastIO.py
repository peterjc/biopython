# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Bio.SearchIO support for NCBI BLAST output."""

from _objects import SearchResult

def BlastXmlIterator(handle):
    from Bio.Blast.NCBIXML import parse
    for index, record in enumerate(parse(handle)):
        #Note this heuristic may need updating for NCBI BLAST 2.2.25+
        query_id = record.query_id
        if str(index+1) == query_id:
            query_id = record.query.split(None,1)[0]
        yield SearchResult(query_id, len(record.alignments))

