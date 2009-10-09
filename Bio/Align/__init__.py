"""Code for dealing with sequence alignments.
"""

from Bio.Align.Generic import Alignment as _BaseAlignment
class MultiSeqAlignment(_BaseAlignment) :
    """A multiple sequence alignment (all entries the same length).

    This class is used to represent a classical multiple sequence alignment,
    for example from ClustalW, PHYLIP or MUSCLE. The key point is that all the
    sequences within it are the same length (but will often include gaps, and
    leading or trailing padding). These objects are used in Bio.AlignIO, and
    in many ways are like a list of SeqRecord objects.

    This class is not really suitable for second generation sequencing contig
    output, which can be regarded as a consensus sequence and many short reads.
    """
    pass
