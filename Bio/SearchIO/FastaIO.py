 # Copyright 2008-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Bio.SearchIO support for "fasta-m10" output from Bill Pearson's FASTA tools.

You are expected to use this module via the Bio.SearchIO or Bio.AlignIO
functions (or the Bio.SeqIO functions if you want to work directly with the
gapped sequences).

This module contains a parser for the pairwise alignments produced by Bill
Pearson's FASTA tools, for use from the Bio.SearchIO and Bio.AlignIO interface
where it is refered to as the "fasta-m10" file format (as we only support the
machine readable output format selected with the -m 10 command line option).

This module does NOT cover the generic "fasta" file format originally
developed as an input format to the FASTA tools.  The Bio.AlignIO and
Bio.SeqIO both use the Bio.SeqIO.FastaIO module to deal with these files,
which can also be used to store a multiple sequence alignments.

Here is an example using Bio.SearchIO to loop over a five query file where
some queries have no hits:

    >>> from Bio import SearchIO
    >>> for result in SearchIO.parse("Fasta/output003.m10", "fasta-m10"):
    ...     print "Query %s has %i matches" % (result.query_id, len(result))
    Query gi|152973837|ref|YP_001338874.1| has 1 matches
    Query gi|152973838|ref|YP_001338875.1| has 0 matches
    Query gi|152973839|ref|YP_001338876.1| has 0 matches
    Query gi|152973840|ref|YP_001338877.1| has 1 matches
    Query gi|152973841|ref|YP_001338878.1| has 1 matches

If we use Bio.AlignIO we don't see the queries with no matches, and since there
is just one HSP per match, we get three alignments only:

    >>> from Bio import AlignIO
    >>> for i, result in enumerate(AlignIO.parse("Fasta/output003.m10", "fasta-m10")):
    ...     print i
    ...     print result
    0
    SingleLetterAlphabet() alignment with 2 rows and 55 columns
    ISISNNKDQYEELQKEQGERDLKTVDQLVRIAAAGGGLRLSAST...IAA gi|152973837|ref|YP_001338874.1|
    VRLTAEEDQ--EIRKRAAECG-KTVSGFLRAAALGKKVNSLTDD...LGA gi|10955263|ref|NP_052604.1|
    1
    SingleLetterAlphabet() alignment with 2 rows and 22 columns
    DDAEHLFRTLSSR-LDALQDGN gi|152973840|ref|YP_001338877.1|
    DDRANLFEFLSEEGITITEDNN gi|10955265|ref|NP_052606.1|
    2
    SingleLetterAlphabet() alignment with 2 rows and 63 columns
    VFGSFEQPKGEHLSGQVSEQ--RDTAFADQNEQVIRHLKQEIEH...QAM gi|152973841|ref|YP_001338878.1|
    VYTSFN---GEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDK...KGI gi|10955264|ref|NP_052605.1|

And here is what it looks like with the low level iterator (which you probably
won't ever need to use):

    >>> from Bio.SearchIO.FastaIO import FastaM10Iterator
    >>> handle = open("Fasta/output003.m10")
    >>> for i, result in enumerate(FastaM10Iterator(handle)):
    ...     print i
    ...     if isinstance(result, basestring):
    ...         print "%s has no matches" % result
    ...     else:
    ...         print result
    0
    SingleLetterAlphabet() alignment with 2 rows and 55 columns
    ISISNNKDQYEELQKEQGERDLKTVDQLVRIAAAGGGLRLSAST...IAA gi|152973837|ref|YP_001338874.1|
    VRLTAEEDQ--EIRKRAAECG-KTVSGFLRAAALGKKVNSLTDD...LGA gi|10955263|ref|NP_052604.1|
    1
    gi|152973838|ref|YP_001338875.1| has no matches
    2
    gi|152973839|ref|YP_001338876.1| has no matches
    3
    SingleLetterAlphabet() alignment with 2 rows and 22 columns
    DDAEHLFRTLSSR-LDALQDGN gi|152973840|ref|YP_001338877.1|
    DDRANLFEFLSEEGITITEDNN gi|10955265|ref|NP_052606.1|
    4
    SingleLetterAlphabet() alignment with 2 rows and 63 columns
    VFGSFEQPKGEHLSGQVSEQ--RDTAFADQNEQVIRHLKQEIEH...QAM gi|152973841|ref|YP_001338878.1|
    VYTSFN---GEKFSSYTLNKVTKTDEYNDLSELSASFFKKNFDK...KGI gi|10955264|ref|NP_052605.1|
    >>> handle.close()

"""

#TODO - Once this is fully working, it can be called by Bio.AlignIO.FastaIO
#whose parser currently ignores queries with no hits.

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO.Interfaces import AlignmentIterator
from Bio.Alphabet import single_letter_alphabet, generic_dna, generic_protein
from Bio.Alphabet import Gapped

from _objects import HSPAlignment

def _extract_alignment_region(alignment_seq_with_flanking, annotation):
    """Helper function for the main parsing code (PRIVATE).

    To get the actual pairwise alignment sequences, we must first
    translate the un-gapped sequence based coordinates into positions
    in the gapped sequence (which may have a flanking region shown
    using leading - characters).  To date, I have never seen any
    trailing flanking region shown in the m10 file, but the
    following code should also cope with that.

    Note that this code seems to work fine even when the "sq_offset"
    entries are prsent as a result of using the -X command line option.
    """
    align_stripped = alignment_seq_with_flanking.strip("-")
    display_start = int(annotation['al_display_start'])
    if int(annotation['al_start']) <= int(annotation['al_stop']):
        start = int(annotation['al_start']) \
              - display_start
        end   = int(annotation['al_stop']) \
              - display_start + 1
    else:
        #FASTA has flipped this sequence...
        start = display_start \
              - int(annotation['al_start'])
        end   = display_start \
              - int(annotation['al_stop']) + 1
    end += align_stripped.count("-")
    assert 0 <= start and start < end and end <= len(align_stripped), \
           "Problem with sequence start/stop,\n%s[%i:%i]\n%s" \
           % (alignment_seq_with_flanking, start, end, annotation)
    return align_stripped[start:end]

def FastaM10Iterator(handle):
    """Alignment iterator for the FASTA tool's pairwise alignment output.

    This is for reading the pairwise alignments output by Bill Pearson's
    FASTA program when called with the -m 10 command line option for machine
    readable output.  For more details about the FASTA tools, see the website
    http://fasta.bioch.virginia.edu/ and the paper:

         W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448

    This class is intended to be used via the Bio.SearchIO and Bio.AlignIO
    functions by specifying the format as "fasta-m10".

    Note that this is not a full blown parser for all the information
    in the FASTA output - for example, most of the header and all of the
    footer is ignored.

    Also note that there can be up to about 30 letters of flanking region
    included in the raw FASTA output as contextual information.  This is NOT
    part of the alignment itself, and is not included in the resulting
    MultipleSeqAlignment objects returned.
    """
    state_PREAMBLE = -1
    state_NONE = 0
    state_QUERY_HEADER = 1
    state_ALIGN_HEADER = 2
    state_ALIGN_QUERY = 3
    state_ALIGN_MATCH = 4
    state_ALIGN_CONS = 5

    def build_hsp():
        assert query_tags, query_tags
        assert match_tags, match_tags
        evalue = align_tags.get("fa_expect", None)
        q = "?" #Just for printing len(q) in debug below
        m = "?" #Just for printing len(m) in debug below
        tool = global_tags.get("tool", "").upper()
        try:
            q = _extract_alignment_region(query_seq, query_tags)
            if tool in ["TFASTX"] and len(match_seq) == len(q):
                m = match_seq
                #Quick hack until I can work out how -, * and / characters
                #and the apparent mix of aa and bp coordindates works.
            else:
                m = _extract_alignment_region(match_seq, match_tags)
            assert len(q) == len(m)
        except AssertionError, err:
            print "Darn... amino acids vs nucleotide coordinates?"
            print tool
            print query_seq
            print query_tags
            print q, len(q)
            print match_seq
            print match_tags
            print m, len(m)
            print handle.name
            raise err
        return HSPAlignment(query_id, match_id, evalue, q, m)

    state = state_PREAMBLE
    query_id = None
    match_id = None
    global_tags = {}
    header_tags = {}
    align_tags = {}
    query_tags = {}
    match_tags = {}
    query_seq = ""
    match_seq = ""
    cons_seq = ""
    for line in handle:
        if ">>>" in line and not line.startswith(">>>"):
            if query_id and match_id:
                #This happens on old FASTA output which lacked an end of
                #query >>><<< marker line.
                yield build_hsp()
            state = state_NONE
            query_id = line[line.find(">>>")+3:].split(None,1)[0]
            match_id = None
            header_tags = {}
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
        elif line.startswith("!! No "):
            #e.g.
            #!! No library sequences with E() < 0.5
            #or on more recent versions,
            #No sequences with E() < 0.05
            assert state == state_NONE
            assert not header_tags
            assert not align_tags
            assert not match_tags
            assert not query_tags
            assert match_id is None
            assert not query_seq
            assert not match_seq
            assert not cons_seq
            if query_id is not None:
                #Yield a string to indicate no hits for this query
                yield query_id
            query_id = None
        elif line.strip() in [">>><<<", ">>>///"]:
            #End of query, possible end of all queries
            if query_id and match_id:
                yield build_hsp()
            state = state_NONE
            query_id = None
            match_id = None
            header_tags = {}
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
        elif line.startswith(">>>"):
            #Should be start of a match!
            assert query_id is not None
            assert line[3:].split(", ",1)[0] == query_id, line
            assert match_id is None
            assert not header_tags
            assert not align_tags
            assert not query_tags
            assert not match_tags
            assert not match_seq
            assert not query_seq
            assert not cons_seq
            state = state_QUERY_HEADER
        elif line.startswith(">>"):
            #Should now be at start of a match alignment!
            if query_id and match_id:
                yield build_hsp()
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
            try:
                match_id, match_descr = line[2:].split(None,1)
            except ValueError:
                match_id = line[2:].strip()
                match_descr = ""
            state = state_ALIGN_HEADER
        elif line.startswith(">--"):
            #End of one HSP
            assert query_id and match_id, line
            yield build_hsp()
            #Clean up read for next HSP
            #but reuse header_tags
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
            state = state_ALIGN_HEADER
        elif line.startswith(">"):
            if state == state_ALIGN_HEADER:
                #Should be start of query alignment seq...
                assert query_id is not None, line
                assert match_id is not None, line
                assert query_id.startswith(line[1:].split(None,1)[0]), line
                state = state_ALIGN_QUERY
            elif state == state_ALIGN_QUERY:
                #Should be start of match alignment seq
                assert query_id is not None, line
                assert match_id is not None, line
                assert match_id.startswith(line[1:].split(None,1)[0]), line
                state = state_ALIGN_MATCH
            elif state == state_NONE:
                #Can get > as the last line of a histogram
                pass
            else:
                assert False, "state %i got %r" % (state, line)
        elif line.startswith("; al_cons"):
            assert state == state_ALIGN_MATCH, line
            state = state_ALIGN_CONS
            #Next line(s) should be consensus seq...
        elif line.startswith("; "):
            key, value = [s.strip() for s in line[2:].split(": ",1)]
            if state == state_QUERY_HEADER:
                header_tags[key] = value
            elif state == state_ALIGN_HEADER:
                align_tags[key] = value
            elif state == state_ALIGN_QUERY:
                query_tags[key] = value
            elif state == state_ALIGN_MATCH:
                match_tags[key] = value
            else:
                assert False, "Unexpected state %r, %r" % (state, line)
        elif state == state_ALIGN_QUERY:
            query_seq += line.strip()
        elif state == state_ALIGN_MATCH:
            match_seq += line.strip()
        elif state == state_ALIGN_CONS:
            cons_seq += line.strip("\n")
        elif state == state_PREAMBLE:
            if line.startswith("#"):
                global_tags["command"] = line[1:].strip()
            elif line.startswith(" version "):
                global_tags["version"] = line[9:].strip()
            elif " compares a " in line:
                global_tags["tool"] = line[:line.find(" compares a ")].strip()
            elif " searches a " in line:
                global_tags["tool"] = line[:line.find(" searches a ")].strip()
        else:
            pass


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
