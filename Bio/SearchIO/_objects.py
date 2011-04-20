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
    def __init__(self, query_id, hits):
        self.query_id = query_id
        for h in hits:
            if not isinstance(h, TopMatches):
                raise TypeError
        self._hits = hits

    def __len__(self):
        return len(self._hits)

    def __iter__(self):
        return iter(self._hits)

class TopMatches(object):
    """Placeholder object to store HSPs/alignments for one match.

    Contains one or more pairwise HSPs/alignments (although these may just
    have co-ordinates and a score, but not the sequences).
    """
    def __init__(self, match_id, alignments):
        self.match_id = match_id
        self._aligns = alignments
        if not alignments:
            raise ValueError("Expect at least one HSP/alignment")

    def __str__(self):
        #return "%s [%i alignments]" % (self.match_id, len(self))
        return self.match_id

    def __len__(self):
        return len(self._aligns)

    def __iter__(self):
        return iter(self._aligns)
