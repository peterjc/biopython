# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Bio.SearchIO objects (private)."""

class SearchResult(object):
    """Placeholder object to store search results for one query.

    Contains zero or more TopMatches objects.
    """
    def __init__(self, query_id, matches):
        self.query_id = query_id
        for match in matches:
            if not isinstance(match, TopMatches):
                raise TypeError("Expected a TopMatches, not %r" % match)
            if match.query_id != query_id:
                raise ValueError("Query ID doesn't agree")
        self._matches = matches

    def __len__(self):
        return len(self._matches)

    def __iter__(self):
        return iter(self._matches)


class TopMatches(object):
    """Placeholder object to store HSPs/alignments for one match.

    Contains one or more pairwise HSPs/alignments (although these may just
    have co-ordinates and a score, but not the sequences).
    """
    def __init__(self, query_id, match_id, hits):
        self.query_id = query_id
        self.match_id = match_id
        self._hits = hits
        if not hits:
            raise ValueError("Expect at least one HSP/alignment")
        for hit in hits:
            if not isinstance(hit, MatchHit):
                raise TypeError("Expected a MatchHit, not %r" % hit)
            if hit.query_id != query_id:
                raise ValueError("Query ID doesn't agree")
            if hit.match_id != match_id:
                raise ValueError("Match ID doesn't agree")

    def __str__(self):
        #return "%s [%i alignments]" % (self.match_id, len(self))
        return self.match_id

    def __len__(self):
        return len(self._hits)

    def __iter__(self):
        return iter(self._hits)

class MatchHit(object):
    """Placeholder object to store a single HSPs/alignment.

    Initially we'll just store the ID and e-value, but a subclass of this
    and the multiple sequence alignment class seems useful for holding a
    pairwise alignment.
    """
    def __init__(self, query_id, match_id, evalue):
        self.query_id = query_id
        self.match_id = match_id
        self.evalue = evalue
