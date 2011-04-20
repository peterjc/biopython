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
        query_id = record.query_id
        if query_id in ["%i" % (index+1), "Query_%i" % (index+1)]:
            query_id = record.query.split(None,1)[0]
        matches = []
        for alignment in record.alignments:
            match_id = alignment.hit_id
            #Expect to need heuristic here too...
            matches.append(match_id)
        yield SearchResult(query_id, matches)

def BlastPairwiseTextIterator(handle):
    from Bio.Blast.NCBIStandalone import BlastParser, Iterator
    for record in Iterator(handle, BlastParser()):
        query_id = record.query_id
        if query_id is None:
            query_id = record.query.split(None,1)[0]
        #Get matches in the table of one line descriptions
        matches_d = [description.title.split(None,1)[0] \
                     for description in record.descriptions]
        #Get matches in the pairwise alignments section
        matches_a = []
        for alignment in record.alignments:
            if alignment.hit_id:
                matches_a.append(alignment.hit_id)
            elif alignment.hit_def:
                matches_a.append(alignment.hit_def)
            elif alignment.title.startswith(">"):
                matches_a.append(alignment.title[1:].split(None,1)[0])
            else:
                matches_a.append(alignment.title.split(None,1)[0])
        #Check they agree
        m = min(len(matches_a), len(matches_d))
        assert matches_a[:m] == matches_d[:m], (query_id, matches_a, matches_d)
        #Take whichever is longer (controlled by different limits)
        if len(matches_a) > len(matches_d):
            matches = matches_a
        else:
            matches = matches_d
        yield SearchResult(query_id, matches)
                            
