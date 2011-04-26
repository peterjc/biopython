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
information. Observe this small BLASTP search using the XML output:

    >>> from Bio import SearchIO
    >>> filename = "Blast/blastp_four_human_vs_rhodopsin.xml"
    >>> for search in SearchIO.parse(filename, "blast-xml"):
    ...     print search.query_id, len(search)
    sp|Q9BS26|ERP44_HUMAN 0
    sp|Q9NSY1|BMP2K_HUMAN 0
    sp|P06213|INSR_HUMAN 0
    sp|P08100|OPSD_HUMAN 6

Now using the same search results but as tabular output:

    >>> from Bio import SearchIO
    >>> filename = "Blast/blastp_four_human_vs_rhodopsin.tabular"
    >>> for search in SearchIO.parse(filename, "blast-stdtab"):
    ...     print search.query_id, len(search)
    sp|P08100|OPSD_HUMAN 6

Notice the tabular output does not contain information about queries with no
matches. In fact, this used to happen with the XML output on older versions
of BLAST too.

Also, the standard 12 column tabular output does not contain the pairwise
alignments of the HSPs either.
"""

from _objects import SearchResult, TopMatches, HSP, HSPAlignment

def _evil_hack(iterator):
    """Copes with BLAST+ XML where querys are wrongly split up (PRIVATE)."""
    prev = None
    for rec in iterator:
        if prev is not None and prev.query_id == rec.query_id:
            prev.alignments.extend(rec.alignments)
            prev.descriptions.extend(rec.descriptions)
        else:
            if prev is not None:
                yield prev
            prev = rec
    if prev is not None:
        yield prev

def BlastXmlIterator(handle):
    from Bio.Blast.NCBIXML import parse
    for index, record in enumerate(_evil_hack(parse(handle))):
        query_id = record.query_id
        if query_id in ["%i" % (index+1), "Query_%i" % (index+1)]:
            query_id = record.query.split(None,1)[0]
        matches = []
        for alignment in record.alignments:
            match_id = alignment.hit_id
            #Expect to need heuristic here too...
            #TODO - Pick sensible alphabet
            hsps = [HSPAlignment(query_id, match_id, h.expect,
                                 h.query, h.sbjct) \
                    for h in alignment.hsps]
            matches.append(TopMatches(query_id, match_id, hsps))
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
            hsps = [HSPAlignment(query_id, match_id, h.expect,
                                 h.query, h.sbjct) \
                    for h in alignment.hsps]
            matches.append(TopMatches(query_id, match_id, hsps))
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
            hsp = HSP(query_id, match_id, parts[10])
            if match_id in [m.match_id for m in matches]:
                matches[-1].append(hsp)
            else:
                matches.append(TopMatches(query_id, match_id, [hsp]))
        else:
            if query_id is not None:
                yield SearchResult(query_id, matches)
            query_id = parts[0]
            hsp = HSP(query_id, match_id, parts[10])
            matches = [TopMatches(query_id, match_id, [hsp])]
    if query_id is not None:
        yield SearchResult(query_id, matches)

def _test():
    """Run the Bio.SearchIO.BlastIO module's doctests (PRIVATE).

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

