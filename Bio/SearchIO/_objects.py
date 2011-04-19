# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Bio.SearchIO objects (private)."""

class SearchResult(object):
    """Placeholder object to store search results for one query."""
    def __init__(self, query_id, hits):
        self.query_id = query_id
        self.hits = hits
    def __len__(self):
        return self.hits #Will be a list later...

