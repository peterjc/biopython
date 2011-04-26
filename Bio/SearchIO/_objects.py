# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Bio.SearchIO objects (private)."""

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

class SearchResult(object):
    """Placeholder object to store search results for one query.

    Contains zero or more TopMatches objects.
    
    Acts a bit liks a list and a dictionary:

    >>> from Bio import SearchIO
    >>> search = SearchIO.read("Blast/xbt006.xml", "blast-xml")
    >>> print len(search)
    10

    Iteration gives you the matches (like a list):

    >>> for match in search:
    ...     print match.match_id, len(match)
    gi|151942244|gb|EDN60600.1| 1
    gi|476059|emb|CAA55606.1| 1
    gi|6320473|ref|NP_010553.1| 1
    gi|61679798|pdb|1R5M|A 1
    gi|6319579|ref|NP_009661.1| 1
    gi|151946495|gb|EDN64717.1| 1
    gi|151943708|gb|EDN62018.1| 2
    gi|151567866|pdb|2PM7|B 2
    gi|6321338|ref|NP_011415.1| 2
    gi|151567870|pdb|2PM9|B 2

    You can access the matches by index (like a list):

    >>> print search[-1].match_id
    gi|151567870|pdb|2PM9|B

    But you can also access the matches by their ID (like a dict):

    >>> print len(search["gi|151567870|pdb|2PM9|B"])
    2

    This is particularly powerful in combination with Bio.SearchIO.index(...)
    to give random access to any matches for any query:

    >>> matches = SearchIO.index("Blast/xbt010.xml", "blast-xml")
    >>> query_id = "gi|4959044|gb|AAD34209.1|AF069992_1"
    >>> match_id = "gi|47078289|ref|NP_035406.3|"
    >>> for hit in matches[query_id][match_id]:
    ...     print hit
    SingleLetterAlphabet() alignment with 2 rows and 600 columns
    MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYR...SVV gi|4959044|gb|AAD34209.1|AF069992_1
    MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYR...SVV gi|47078289|ref|NP_035406.3|

    If a query had no hits against a particular match, you'd get a KeyError
    (like a dict).
    """
    def __init__(self, query_id, matches):
        self.query_id = query_id
        for match in matches:
            if not isinstance(match, TopMatches):
                raise TypeError("Expected a TopMatches, not %r" % match)
            if match.query_id != query_id:
                raise ValueError("Query ID doesn't agree")
        self._matches = matches

    def __repr__(self):
        return "SearchResult(%r, %i matches)" % (self.query_id, len(self))

    def __len__(self):
        return len(self._matches)

    def __iter__(self):
        return iter(self._matches)
    
    def __getitem__(self, value):
        if isinstance(value, int):
            return self._matches[value]
        else:
            for m in self:
                if m.match_id == value:
                    return m
            raise KeyError(value)
    
    def append(self, item):
        """Add a TopMatches or HSP object to the results."""
        try:
            if item.query_id != self.query_id:
                raise ValueError("Query ID %r does not equal %r" \
                                 % (item.query_id, self.query_id))
        except AttributeError:
            raise TypeError("Object should have query_id, e.g. HSP")
        if isinstance(item, TopMatches):
            for m in self:
                if item.match_id == m.match_id:
                    raise ValueError("Match ID %r already present" \
                                     % item.match_id)
            self._matches.append(item)
        elif isinstance(item, HSP):
            for m in self:
                if item.match_id == m.match_id:
                    m.append(item)
                    return
            #Must create container...
            self._matches.append(TopMatches(item.query_id,
                                            item.match_id,
                                            [item]))
        else:
            raise TypeError
    

class TopMatches(object):
    """Placeholder object to store HSPs (pairwise alignments) for one match.

    Contains one or more pairwise HSPs/alignments (although these may just
    have co-ordinates and a score, but not the sequences).
    """
    def __init__(self, query_id, match_id, hsps):
        self.query_id = query_id
        self.match_id = match_id
        self._hsps = hsps
        if not hsps:
            raise ValueError("Expect at least one HSP (pairwise alignment)")
        for hsp in hsps:
            if not isinstance(hsp, HSP):
                raise TypeError("Expected an HSP, not %r" % hit)
            if hsp.query_id != query_id:
                raise ValueError("Query ID doesn't agree")
            if hsp.match_id != match_id:
                raise ValueError("Match ID doesn't agree")

    def __str__(self):
        #return "%s [%i alignments]" % (self.match_id, len(self))
        return self.match_id

    def __repr__(self):
        return "TopMatches(%r, %r, %i hits)" \
               % (self.query_id, self.match_id, len(self))

    def __len__(self):
        return len(self._hsps)

    def __iter__(self):
        return iter(self._hsps)

    def append(self, item):
        """Add an HSP object to the results."""
        if not isinstance(item, HSP):
            raise TypeError
        if item.query_id != self.query_id:
            raise ValueError("Query ID %r does not equal %r" \
                             % (item.query_id, self.query_id))
        if item.match_id != self.match_id:
            raise ValueError("Match ID %r does not equal %r" \
                             % (item.match_id, self.match_id))
        self._hsps.append(item)
 

class HSP(object):
    """High-scoring segment pair (HSP), without sequence.

    Initially we'll just store the ID and e-value, but there is a subclass of
    this and the MultipleSeqAlignment class for holding a pairwise alignment.
    """
    def __init__(self, query_id, match_id, evalue):
        self.query_id = query_id
        self.match_id = match_id
        self.evalue = evalue

    def __repr__(self):
        return "MatchHit(%r, %r, %r)" \
               % (self.query_id, self.match_id, self.evalue)

class HSPAlignment(MultipleSeqAlignment, HSP):
    """High-scoring segment pair (HSP), with aligned sequences.
    
    This is both an HSP object from Bio.SearchIO, and a pairwise
    MultipleSeqAlignment object which can be used with Bio.AlignIO too.
    """
    #TODO - What about general annotation on a MultipleSeqAlignment?
    def __init__(self, query_id, match_id, evalue,
                 aligned_query, aligned_match,
                 alphabet=single_letter_alphabet):
        HSP.__init__(self, query_id, match_id, evalue)
        query = SeqRecord(Seq(aligned_query, alphabet),
                          id=query_id, name="query",
                          description="")
        match = SeqRecord(Seq(aligned_match, alphabet),
                          id=match_id, name="match",
                          description="")
        MultipleSeqAlignment.__init__(self, [query, match], alphabet)

def _test():
    """Run the Bio.SearchIO._object module's doctests (PRIVATE).

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..","..","Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","..","Tests"))
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
