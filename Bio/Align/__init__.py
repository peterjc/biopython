"""Code for dealing with sequence alignments.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Alphabet
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
    def __init__(self, records, alphabet=Alphabet.single_letter_alphabet) :
        #Temporary, going to allow a list of records to be given...
        if isinstance(records, Alphabet.Alphabet) \
        or isinstance(records, Alphabet.AlphabetEncoder):
            records, alphabet = [], records
            import warnings
            warnings.warn("Bio.Align.MultiSeqAlignment takes different "
                          "arguments to the older Bio.Align.Generic.Alphabet.",
                          DeprecationWarning)
        if not (isinstance(alphabet, Alphabet.Alphabet) \
        or isinstance(alphabet, Alphabet.AlphabetEncoder)):
            raise ValueError("Invalid alphabet argument")
        self._alphabet = alphabet
        # hold everything at a list of SeqRecord objects
        self._records = []
        self._length = 0
        self.extend(records)
        #TODO - If the provided alphabet was the default, should
        #we pick a consensus based on the records (if any)?

    def _check_len(self) :
        """Debugging method (PRIVATE)."""
        for record in self._records :
            assert self._length == len(record), \
                   "self._length = %i but record length %i\n%s" \
                   % (self._length, len(record), str(self._records))
        
    def get_alignment_length(self) :
        """Returns the sequence length of the records in the alignment."""
        if not self._length and self._records and len(self._records[0]):
            #Can happen if the private list self._records is edited directly!
            raise RuntimeError("Private variables _records and _length out of sync!")
        self._check_len()
        return self._length

    def append(self, record) :
        """Add another SeqRecord object as the last row of the alignment."""
        if not isinstance(record, SeqRecord) :
            raise TypeError
        if self._length and len(record) != self._length :
            raise ValueError("SeqRecord must be length %i" % self._length)
        #TODO - Reuse Bio.SeqIO.to_alignment() code to check the alphabet
        self._records.append(record)
        if not self._length :
            self._length = len(record)
        
    def extend(self, records) :
        """Add more SeqRecord objects as row of the alignment."""
        for record in records :
            self.append(record)
        self._check_len()

    def add_sequence(self, descriptor, sequence, start = None, end = None,
                     weight = 1.0):
        if start or end or weight != 1.0 :
            raise ValueError("This method exists only to provide limited "
                      "backwards compatibility with the older "
                      "Bio.Align.Generic.Alphabet class. However, it does not "
                      "support the start, end and weight parameters.")
        self.append(SeqRecord(seq=Seq(sequence, self._alphabet),
                              id=descriptor,
                              description=descriptor))
