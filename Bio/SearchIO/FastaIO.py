# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Bio.SearchIO support for "fasta-m10" output from Bill Pearson's FASTA tools."""

#TODO - Once this is fully working, it can be called by Bio.AlignIO.FastaIO
#whose parser currently ignores queries with no hits.

from _objects import SearchResult, TopMatches, HSP


def FastaM10Iterator(handle):
    line = handle.readline()
    while line and ">>><<<" not in line:
        while line.startswith("#") or ">>>" not in line:
            #Skip the file header before the alignments.
            line = handle.readline()
        assert ">>>" in line and not line.startswith(">>>"), line
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
        #Sometimes have queries with no matches:
        """
          2>>>gi|152973838|ref|YP_001338875.1| hypothetical protein KPN_pKPN7p10263 [K... pneumonia 76 aa - 76 aa
         vs  NC_002127.faa library

            579 residues in     3 sequences
         Altschul/Gish params: n0: 76 Lambda: 0.158 K: 0.019 H: 0.100

        FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
         join: 36, opt: 24, open/ext: -10/-2, width:  16
         Scan time:  0.000
        !! No library sequences with E() < 0.5
        """
        query_id = line.split(">>>",1)[1].split(None,1)[0]
        matches = []
        while True:
            line = handle.readline()
            if line.startswith("!! No library sequences with"):
                assert not matches, (query_id, matches)
                break
            elif ">>><<<" in line or not line:
                #End of file!
                break
            elif line.startswith(">>>"):
                assert line.startswith(">>>"+query_id), line
            elif ">>>" in line:
                #Next query!
                break
            elif line.startswith(">>"):
                match_id = line[2:].split(None,1)[0]
                hsps = [HSP(query_id, match_id, "e?")]
                matches.append(TopMatches(query_id, match_id, hsps))
        yield SearchResult(query_id, matches)

