from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class LazySeqRecord(SeqRecord):

    def __init__(self, handle, offset, raw_len,
                 id, index, alphabet, seq_len=None):
        """Create LazySeqRecord which will access file on demand.

        handle - input file handle
        offset - position in the handle where this record starts
        raw_len - length of raw record in the file (optional?)
        id - string to be used as the SeqRecord's id
        index - slice object describing which part of sequence wanted
        alphabet - alphabet to use for the sequence (optional?)
        seq_len - length of the (full) sequence, see index argument (optional)

        The idea of having the index as part of the LazySeqRecord is to allow
        fast slicing and (in principle) this gives us a way to be clever about
        slicing the sequence, per-letter-annotation or features. e.g. When
        asked to load the features, we only need to return those within the
        slice. This is aimed to fit in with indexed access to big files, such
        as SAM/BAM (see samtools), GFF (perhaps using tabix, see samtools), or
        even FASTA in principle (see samtools FASTA indexing which requires
        uniform line wrapping).
        """
        #Note that the id string is supplied as part of the __init__ args for
        #two reasons. First, this is probably the most accessed attribute, so
        #there is little benefit to making it lazy (hunch). Second, I want to
        #reuse the iteration support developed in Bio/SeqIO/_index.py which
        #already returns id, offset, raw_len tuples.
        #
        #NOTE: The subclass implementations of __init__ should not access the
        #handle (for speed). In practice if using a for loop (etc) with a lazy
        #SeqIO parser, the user will access methods of the current record which
        #will move the read position. This means the iterator has to cope with
        #the read position moving during the each yield statement. We could use
        #a separate handle for the lazy records and the parser I suppose... but
        #one handle for each lazy record is a bad idea.
        #
        #TODO - Computing the sequence length up front is often expensive
        #(slow), but if access this on demand needs more infrastructure
        #to handle start/end for __getitem__, __len__, etc.
        self._handle = handle
        self._offset = offset
        self._raw_len = raw_len
        self.id = id
        self._seq_len = seq_len
        if index.step != 1 and index.step is not None:
            raise ValueError
        self._index = index
        self._alphabet = alphabet

    def toseqrecord(self):
        """Turn the LazySeqRecord into a traditional in memory SeqRecord."""
        return SeqRecord(self.seq, self.id, self.name, self.description,
                         self.dbxrefs, self.features, self.annotations,
                         self.letter_annotations)
    
    def __len__(self):
        if self._seq_len is None:
            self._seq_len = self._load_full_len()
        start, end, step = self._index.indices(self._seq_len)
        assert step == 1
        return (end - start) // step

    def __getitem__(self, index):
        if isinstance(index, slice):
            if self._seq_len is None:
                self._seq_len = self._load_len()
            start, stop, step = self._index.indices(self._seq_len)
            if index.step is None or index.step == 1:
                if index.start >= 0:
                    new_start = min(start + index.start, self._seq_len)
                else:
                    new_start = max(0, self._stop + index.start)
                if index.stop >= 0:
                    new_stop = min(start + index.stop, self._seq_len)
                else:
                    new_stop = max(0, stop + index.stop)
                return self.__class__(self._handle,
                                      self._offset, self._raw_len,
                                      self.id, slice(new_start, new_stop),
                                      self._alphabet, self._seq_len)
            else:
                #Can we cope with a step as well? Problem is slice of slice...
                return self.toseqrecord()[index]
        elif isinstance(index, int):
            #Extract single letter of sequence...
            return self.seq[index]
        else:
            raise ValueError, "Invalid index"
    
    def _get_seq(self):
        try:
            return self._seq
        except AttributeError:
            #Load it now
            s = Seq(self._load_seq(self._index), self._alphabet)
            self._seq = s
            return s
    seq = property(fget=_get_seq,
                   fset=SeqRecord._set_seq, #Defined in base classes
                   doc="The sequence itself, as a Seq or MutableSeq object.")

    def _get_name(self):
        try:
            return self._name
        except AttributeError:
            #Load it now
            temp = self._load_name()
            self._id = temp
            return temp
    def _set_name(self, value):
        self._name = value
    name = property(fget=_get_name,
                    fset=_set_name,
                    doc="The sequence name (string)")

    def _get_description(self):
        try:
            return self._description
        except AttributeError:
            #Load it now
            temp = self._load_description()
            self._description = temp
            return temp
    def _set_description(self, value):
        self._description = value
    description = property(fget=_get_description,
                           fset=_set_description,
                           doc="The sequence description (string)")

    def _get_annotations(self):
        try:
            return self._annotations
        except AttributeError:
            #Load it now
            temp = self._load_annotations()
            self._annotations = temp
            return temp
    def _set_annotations(self, value):
        if not isinstance(value, dict):
            raise TypeError
        self._annotations = value
    annotations = property(fget=_get_annotations,
                           fset=_set_annotations,
                           doc="General annotation (dict)")

    def _get_per_letter_annotations(self):
        try:
            return self._per_letter_annotations
        except AttributeError:
            #Load it now
            temp = self._load_per_letter_annotations(self._index)
            SeqRecord._set_per_letter_annotations(self, temp) #validates it!
            return temp
    letter_annotations = property( \
        fget=_get_per_letter_annotations,
        fset=SeqRecord._set_per_letter_annotations,
        doc = SeqRecord.letter_annotations.__doc__)

    def _get_features(self):
        try:
            return self._features
        except AttributeError:
            #Load it now
            temp = self._load_features(self._index)
            self._features = temp
            return temp
    def _set_features(self, value):
        if not isinstance(value, list):
            raise TypeError
        self._features = value
    features = property(fget=_get_features,
                        fset=_set_features,
                        doc="Features (list of SeqFeature objects)")

    def _get_dbxrefs(self):
        try:
            return self._dbxrefs
        except AttributeError:
            #Load it now
            if self._seq_len is None:
                self._seq_len = self._load_len()
            if (0, self._seq_len, 1) == self._index.indices(self._seq_len):
                temp = self._load_dbxrefs()
            else:
                #Slicing the SeqRecord would discard the cross references,
                temp = []
            self._dbxrefs = temp
            return temp
    def _set_dbxrefs(self, value):
        if not isinstance(value, list):
            raise TypeError
        self._dbxrefs = value
    dbxrefs = property(fget=_get_dbxrefs,
                       fset=_set_dbxrefs,
                       doc="Database cross-references (list)")

    def _load_full_len(self):
        """Extracts the full sequence length from the file."""
        #Default implementation, bit of a hack.
        #For some formats this can be done very efficiently.
        #e.g. GenBank files should have length in their LOCUS line
        return len(self._load_seq(slice(None,None)))
        
    def _load_seq(self, index):
        """Extracts the (sub)sequence from the file."""
        raise NotImplementedError

    def _load_name(self):
        """Extracts name from the file."""
        return ""

    def _load_description(self):
        """Extracts description from the file."""
        return ""

    def _load_annotations(self):
        """Extracts annotations from the file."""
        return {}

    def _load_per_letter_annotations(self, index):
        """Extracts per-letter-annotation for (sub)sequence from the file."""
        return {}

    def _load_features(self, index):
        """Extracts features for (sub)sequence from the file."""
        return []

    def _load_dbxrefs(self):
        """Extracts database cross references from the file."""
        return []
