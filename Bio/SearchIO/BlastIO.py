# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Bio.SearchIO support for NCBI BLAST output.

Rather than accessing the Bio.SearchIO.BlastIO module directly, you are
expected to access its functionality via the top level functions of
Bio.SearchIO under the format names "blast-xml", "blast-stdtab" and
"blast-text".

Note that these different BLAST output formats contain somewhat different
information. In particular, the standard 12 column tabular output does not
contain information about queries with no hits, and does not contain the
pairwise alignments of the hits either.
"""

from _objects import SearchResult, TopMatches, MatchHit, MatchHitAlignment

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
            #TODO - Pick sensible alphabet
            hits = [MatchHitAlignment(query_id, match_id, h.expect,
                                      h.query, h.sbjct) \
                    for h in alignment.hsps]
            matches.append(TopMatches(query_id, match_id, hits))
        yield SearchResult(query_id, matches)

def BlastPairwiseTextIterator(handle):
    from Bio.Blast.NCBIStandalone import BlastParser, Iterator
    for record in Iterator(handle, BlastParser()):
        query_id = record.query_id
        if query_id is None:
            query_id = record.query.split(None,1)[0]
        #Get matches in the pairwise alignments section
        #(Ignore the separate one line descriptions given before this)
        matches = []
        for alignment in record.alignments:
            if alignment.hit_id:
                match_id = alignment.hit_id
            elif alignment.hit_def:
                match_id = alignment.hit_def
            elif alignment.title.startswith(">"):
                match_id = alignment.title[1:].split(None,1)[0]
            else:
                match_id = alignment.title.split(None,1)[0]
            #TODO - Pick sensible alphabet
            hits = [MatchHitAlignment(query_id, match_id, h.expect,
                                      h.query, h.sbjct) \
                    for h in alignment.hsps]
            matches.append(TopMatches(query_id, match_id, hits))
        yield SearchResult(query_id, matches)
                            
def BlastStandardTabularIterator(handle):
    query_id = None
    matches = []
    for line in handle:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            raise ValueError("Only %i columns in line %r" \
                             % (len(parts), line))
        match_id = parts[1]
        if query_id == parts[0]:
            hit = MatchHit(query_id, match_id, parts[10])
            if match_id in [m.match_id for m in matches]:
                #Need a proper way to add an HSP/alignment...
                #...for now we're just storing the e-value
                matches[-1]._hits.append(hit)
            else:
                matches.append(TopMatches(query_id, match_id, [hit]))
        else:
            if query_id is not None:
                yield SearchResult(query_id, matches)
            query_id = parts[0]
            hit = MatchHit(query_id, match_id, parts[10])
            matches = [TopMatches(query_id, match_id, [hit])]
    if query_id is not None:
        yield SearchResult(query_id, matches)
