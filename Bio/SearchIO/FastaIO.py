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


# TODO - Turn this into a doctest
class FastaM10Iterator(AlignmentIterator):
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
    
    def next(self):
        """Reads from the handle to construct and return the next alignment.

        This returns the pairwise alignment of query and match/library
        sequences as an HSPAlignment object (a subclass of the
        MultipleSeqAlignment object) which containings two rows.

        If the next query has no results however, a string is returned (the
        query ID).
        """
        handle = self.handle

        try:
            #Header we saved from when we were parsing
            #the previous alignment.
            line = self._header
            del self._header
        except AttributeError:      
            line = handle.readline()
        if not line:
            raise StopIteration

        if line.startswith("#") or line.startswith(">>><<<"):
            #Skip the file header before the alignments.
            line = self._skip_file_header(line)
        assert not line.startswith(">>><<<")
        while ">>>" in line and not line.startswith(">>>"):
            #Moved onto the next query sequence!
            self._query_descr = ""
            self._query_header_annotation = {}
            #Read in the query header
            line = self._parse_query_header(line)
            #Now should be some alignments, but if not we move onto the next query
            if line.startswith("!! No library sequences with"):
                return self._query_descr.split(None,1)[0] 
        if not line:
            #End of file
            raise StopIteration
        if ">>>///" in line:
            #Reached the end of the alignments, no need to read the footer...
            #Note this is a new line type in fasta-36.3, previously >>><<<
            #was used but that not occurs at the end of each query.
            raise StopIteration
        assert not line.startswith(">>><<<")


        #Should start >>... and not >>>...
        assert line[0:2] == ">>" and not line[2] == ">", repr(line)

        query_seq_parts, match_seq_parts = [], []
        query_annotation, match_annotation = {}, {}
        match_descr = ""
        alignment_annotation = {}

        #This should be followed by the target match ID line, then more tags.
        #e.g.
        """
        >>gi|152973545|ref|YP_001338596.1| putative plasmid SOS inhibition protein A [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
        ; fa_frame: f
        ; fa_initn:  52
        ; fa_init1:  52
        ; fa_opt:  70
        ; fa_z-score: 105.5
        ; fa_bits: 27.5
        ; fa_expect:  0.082
        ; sw_score: 70
        ; sw_ident: 0.279
        ; sw_sim: 0.651
        ; sw_overlap: 43
        """
        if (not line[0:2] == ">>") or line[0:3] == ">>>":
            raise ValueError("Expected target line starting '>>'")
        match_descr = line[2:].strip()
        #Handle the following "alignment hit" tagged data, e.g.
        line = handle.readline()
        line = self._parse_tag_section(line, alignment_annotation)
        assert not line[0:2] == "; "
        
        #Then we have the alignment numbers and sequence for the query
        """
        >gi|10955265| ..
        ; sq_len: 346
        ; sq_offset: 1
        ; sq_type: p
        ; al_start: 197
        ; al_stop: 238
        ; al_display_start: 167
        DFMCSILNMKEIVEQKNKEFNVDIKKETIESELHSKLPKSIDKIHEDIKK
        QLSC-SLIMKKIDVEMEDYSTYCFSALRAIEGFIYQILNDVCNPSSSKNL
        GEYFTENKPKYIIREIHQET
        """
        if not (line[0] == ">" and line.strip().endswith("..")):
            raise ValueError("Expected line starting '>' and ending '..'")
        assert self._query_descr.startswith(line[1:].split(None,1)[0])
        
        #Handle the following "query alignment" tagged data
        line = handle.readline()
        line = self._parse_tag_section(line, query_annotation)
        assert not line[0:2] == "; "

        #Now should have the aligned query sequence (with leading flanking region)
        while not line[0] == ">":
            query_seq_parts.append(line.strip())
            line = handle.readline()
        
        #Handle the following "match alignment" data
        """
        >gi|152973545|ref|YP_001338596.1| ..
        ; sq_len: 242
        ; sq_type: p
        ; al_start: 52
        ; al_stop: 94
        ; al_display_start: 22
        IMTVEEARQRGARLPSMPHVRTFLRLLTGCSRINSDVARRIPGIHRDPKD
        RLSSLKQVEEALDMLISSHGEYCPLPLTMDVQAENFPEVLHTRTVRRLKR
        QDFAFTRKMRREARQVEQSW
        """
        #Match identifier
        if not (line[0] == ">" and line.strip().endswith("..")):
            raise ValueError("Expected line starting '>' and ending '..', got '%s'" % repr(line))
        assert match_descr.startswith(line[1:].split(None,1)[0])
        
        #Tagged data,
        line = handle.readline()
        line = self._parse_tag_section(line, match_annotation)
        assert not line[0:2] == "; "

        #Now should have the aligned query sequence with flanking region...
        #but before that, since FASTA 35.4.1 there can be an consensus here,
        """
        ; al_cons:
        .::. : :. ---.  :: :. . :  ..-:::-:  :.:  ..:...: 
        etc
        """
        while not (line[0:2] == "; " or line[0] == ">" or ">>>" in line):
            match_seq_parts.append(line.strip())
            line = handle.readline()
        if line[0:2] == "; ":
            assert line.strip() == "; al_cons:"
            align_consensus_parts = []
            line = handle.readline()
            while not (line[0:2] == "; " or line[0] == ">" or ">>>" in line):
                align_consensus_parts.append(line.strip())
                line = handle.readline()
            #If we do anything with this in future, must remove any flanking region.
            align_consensus = "".join(align_consensus_parts)
            del align_consensus_parts
            assert not line[0:2] == "; "
        else:
            align_consensus = None
        assert (line[0] == ">" or ">>>" in line)
        self._header = line

        #We built a list of strings and then joined them because
        #its faster than appending to a string.
        query_seq = "".join(query_seq_parts)
        match_seq = "".join(match_seq_parts)
        del query_seq_parts, match_seq_parts
        #Note, query_seq and match_seq will usually be of different lengths, apparently
        #because in the m10 format leading gaps are added but not trailing gaps!

        #Remove the flanking regions,
        query_align_seq = self._extract_alignment_region(query_seq, query_annotation)
        match_align_seq = self._extract_alignment_region(match_seq, match_annotation)
        #How can we do this for the (optional) consensus?

        #The "sq_offset" values can be specified with the -X command line option.
        #They appear to just shift the origin used in the calculation of the coordinates.
        
        if len(query_align_seq) != len(match_align_seq):
            raise ValueError("Problem parsing the alignment sequence coordinates, " 
                             "following should be the same length but are not:\n"
                             "%s - len %i\n%s - len %i" % (query_align_seq,
                                                           len(query_align_seq),
                                                           match_align_seq,
                                                           len(match_align_seq)))
        if "sw_overlap" in alignment_annotation:
            if int(alignment_annotation["sw_overlap"]) != len(query_align_seq):
                raise ValueError("Specified sw_overlap = %s does not match expected value %i" \
                                 % (alignment_annotation["sw_overlap"],
                                    len(query_align_seq)))

        #TODO - Look at the "sq_type" to assign a sensible alphabet?
        alphabet = self.alphabet
        alignment = MultipleSeqAlignment([], alphabet)

        #TODO - Introduce an annotated alignment class?
        #For now, store the annotation a new private property:
        alignment._annotations = {}
        
        #Want to record both the query header tags, and the alignment tags.
        for key, value in self._query_header_annotation.iteritems():
            alignment._annotations[key] = value
        for key, value in alignment_annotation.iteritems():
            alignment._annotations[key] = value
        
        #Query
        #=====
        record = SeqRecord(Seq(query_align_seq, alphabet),
                           id = self._query_descr.split(None,1)[0].strip(","),
                           name = "query",
                           description = self._query_descr,
                           annotations = {"original_length" : int(query_annotation["sq_len"])})
        #TODO - handle start/end coordinates properly. Short term hack for now:
        record._al_start = int(query_annotation["al_start"])
        record._al_stop = int(query_annotation["al_stop"])
        alignment.append(record)

        #TODO - What if a specific alphabet has been requested?
        #TODO - Use an IUPAC alphabet?
        #TODO - Can FASTA output RNA?
        if alphabet == single_letter_alphabet and "sq_type" in query_annotation:
            if query_annotation["sq_type"] == "D":
                record.seq.alphabet = generic_dna
            elif query_annotation["sq_type"] == "p":
                record.seq.alphabet = generic_protein
        if "-" in query_align_seq:
            if not hasattr(record.seq.alphabet,"gap_char"):
                record.seq.alphabet = Gapped(record.seq.alphabet, "-")

        #Match
        #=====
        record = SeqRecord(Seq(match_align_seq, alphabet),
                           id = match_descr.split(None,1)[0].strip(","),
                           name = "match",
                           description = match_descr,
                           annotations = {"original_length" : int(match_annotation["sq_len"])})
        #TODO - handle start/end coordinates properly. Short term hack for now:
        record._al_start = int(match_annotation["al_start"])
        record._al_stop = int(match_annotation["al_stop"])
        alignment.append(record)

        #This is still a very crude way of dealing with the alphabet:
        if alphabet == single_letter_alphabet and "sq_type" in match_annotation:
            if match_annotation["sq_type"] == "D":
                record.seq.alphabet = generic_dna
            elif match_annotation["sq_type"] == "p":
                record.seq.alphabet = generic_protein
        if "-" in match_align_seq:
            if not hasattr(record.seq.alphabet,"gap_char"):
                record.seq.alphabet = Gapped(record.seq.alphabet, "-")

        #return alignment
        query_id = self._query_descr.split(None,1)[0].strip(",")
        match_id = match_descr.split(None,1)[0].strip(",")
        evalue = query_annotation.get("fa_expect", None)
        return HSPAlignment(query_id, match_id, evalue,
                            query_align_seq, match_align_seq,
                            alphabet)

    def _skip_file_header(self, line):
        """Helper function for the main parsing code.

        Skips over the file header region.
        """
        #e.g. This region:
        """
        # /home/xxx/Downloads/FASTA/fasta-35.3.6/fasta35 -Q -H -E 1 -m 10 -X "-5 -5" NC_002127.faa NC_009649.faa
        FASTA searches a protein or DNA sequence data bank
         version 35.03 Feb. 18, 2008
        Please cite:
         W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448

        Query: NC_002127.faa
        """
        #Note that there is no point recording the command line here
        #from the # line, as it is included again in each alignment
        #under the "pg_argv" tag.  Likewise for the program version.
        while line and (">>>" not in line or line.startswith(">>><<<")):
            line = self.handle.readline()
        return line

    def _parse_query_header(self, line):
        """Helper function for the main parsing code.

        Skips over the free format query header, extracting the tagged parameters.

        If there are no hits for the current query, it is skipped entirely."""
        #e.g. this region (where there is often a histogram too):
        """
          2>>>gi|10955264|ref|NP_052605.1| hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa - 126 aa
        Library: NC_009649.faa   45119 residues in   180 sequences

          45119 residues in   180 sequences
        Statistics: (shuffled [500]) Expectation_n fit: rho(ln(x))= 5.0398+/-0.00968; mu= 2.8364+/- 0.508
         mean_var=44.7937+/-10.479, 0's: 0 Z-trim: 0  B-trim: 0 in 0/32
         Lambda= 0.191631
        Algorithm: FASTA (3.5 Sept 2006) [optimized]
        Parameters: BL50 matrix (15:-5) ktup: 2
         join: 36, opt: 24, open/ext: -10/-2, width:  16
         Scan time:  0.040

        The best scores are:                                      opt bits E(180)
        gi|152973462|ref|YP_001338513.1| hypothetical prot ( 101)   58 23.3    0.22
        gi|152973501|ref|YP_001338552.1| hypothetical prot ( 245)   55 22.5    0.93
        """
        #Sometimes have queries with no matches, in which case we continue to the
        #next query block:
        """
          2>>>gi|152973838|ref|YP_001338875.1| hypothetical protein KPN_pKPN7p10263 [Klebsiella pneumoniae subsp. pneumonia 76 aa - 76 aa
         vs  NC_002127.faa library

            579 residues in     3 sequences
         Altschul/Gish params: n0: 76 Lambda: 0.158 K: 0.019 H: 0.100

        FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
         join: 36, opt: 24, open/ext: -10/-2, width:  16
         Scan time:  0.000
        !! No library sequences with E() < 0.5
        """

        self._query_header_annotation = {}
        self._query_descr = ""

        assert ">>>" in line and not line[0:3] == ">>>"
        #There is nothing very useful in this line, the query description is truncated.
        #However, we have to use this for the no hits case.
        self._query_descr = line[line.find(">>>")+3:].split(None,1)[0]
        
        line = self.handle.readline()
        #We ignore the free form text...
        while line[0:3] != ">>>" or line.startswith(">>><<<"):
            #print "Ignoring %s" % line.strip()
            line = self.handle.readline()
            if not line:
                return ""
            if line.startswith("!! No library sequences with"):
                return line

        #Now want to parse this section:
        """
        >>>gi|10955264|ref|NP_052605.1|, 126 aa vs NC_009649.faa library
        ; pg_name: /home/pjcock/Downloads/FASTA/fasta-35.3.6/fasta35
        ; pg_ver: 35.03
        ; pg_argv: /home/pjcock/Downloads/FASTA/fasta-35.3.6/fasta35 -Q -H -E 1 -m 10 -X -5 -5 NC_002127.faa NC_009649.faa
        ; pg_name: FASTA
        ; pg_ver: 3.5 Sept 2006
        ; pg_matrix: BL50 (15:-5)
        ; pg_open-ext: -10 -2
        ; pg_ktup: 2
        ; pg_optcut: 24
        ; pg_cgap: 36
        ; mp_extrap: 60000 500
        ; mp_stats: (shuffled [500]) Expectation_n fit: rho(ln(x))= 5.0398+/-0.00968; mu= 2.8364+/- 0.508  mean_var=44.7937+/-10.479, 0's: 0 Z-trim: 0  B-trim: 0 in 0/32  Lambda= 0.191631
        ; mp_KS: -0.0000 (N=1066338402) at  20
        ; mp_Algorithm: FASTA (3.5 Sept 2006) [optimized]
        ; mp_Parameters: BL50 matrix (15:-5) ktup: 2  join: 36, opt: 24, open/ext: -10/-2, width:  16
        """

        assert not line.startswith(">>><<<"), line
        assert line[0:3] == ">>>", line
        self._query_descr = line[3:].strip()

        #Handle the following "program" tagged data,
        line = self.handle.readline()
        line = self._parse_tag_section(line, self._query_header_annotation)
        assert not line[0:2] == "; ", line
        assert line[0:2] == ">>" or ">>>" in line, line
        return line


    def _extract_alignment_region(self, alignment_seq_with_flanking, annotation):
        """Helper function for the main parsing code.

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
                  - display_start \
                  + align_stripped.count("-") + 1
        else:
            #FASTA has flipped this sequence...
            start = display_start \
                  - int(annotation['al_start'])
            end   = display_start \
                  - int(annotation['al_stop']) \
                  + align_stripped.count("-") + 1
        assert 0 <= start and start < end and end <= len(align_stripped), \
               "Problem with sequence start/stop,\n%s[%i:%i]\n%s" \
               % (alignment_seq_with_flanking, start, end, annotation)
        return align_stripped[start:end]


    def _parse_tag_section(self, line, dictionary):
        """Helper function for the main parsing code.

        line - supply line just read from the handle containing the start of
               the tagged section.
        dictionary - where to record the tagged values

        Returns a string, the first line following the tagged section."""
        if not line[0:2] == "; ":
            raise ValueError("Expected line starting '; '")
        while line[0:2] == "; ":
            tag, value = line[2:].strip().split(": ",1)
            #fasta34 and early versions of fasta35 will reuse the pg_name and
            #pg_ver tags for the program executable and name, and the program
            #version and the algorithm version, respectively.  This is fixed
            #in FASTA 35.4.1, but we can't assume the tags are unique:
            #if tag in dictionary:
            #    raise ValueError("Repeated tag '%s' in section" % tag)
            dictionary[tag] = value
            line = self.handle.readline()
        return line


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
