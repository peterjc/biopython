from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class LazySeqRecord(SeqRecord):
    _raw_formats = [] #optimisation for the format method

    def __init__(self, handle, offset, raw_len, seq_len, index, alphabet):
        self._handle = handle
        self._offset = offset
        self._raw_len = raw_len
        self._seq_len = seq_len
        self._start, self._stop, step = index.indices(seq_len)
        if step != 1:
            raise ValueError
        self._alphabet = alphabet

    def __format__(self, format_spec):
        #Can we just read the raw read from the disk?
        if format_spec in self._raw_formats and self._raw_len:
            print "Speedy %s" % format_spec
            h = self._handle
            h.seek(self._offset)
            return h.read(self._raw_len)
        else:
            print "Slow %s, not in %r" % (format_spec, self._raw_formats)
            return SeqRecord.__format__(self, format_spec)
    
    def __len__(self):
        #start, end, step = self._index.indices(self._seq_len)
        #return (end - step) // step
        l = self._stop - self._start
        assert l >= 0, l
        return l

    def __getitem__(self, index):
        if isinstance(index, slice):
            if index.step is None or index.step == 1:
                if index.start >= 0:
                    new_start = min(self._start + index.start, self._seq_len)
                else:
                    new_start = max(0, self._stop + index.start)
                if index.stop >= 0:
                    new_stop = min(self._start + index.stop, self._seq_len)
                else:
                    new_stop = max(0, self._stop + index.stop)
                return self.__class__(self._handle, self._offset, self._raw_len,
                                      self._seq_len, slice(new_start, new_stop),
                                      self._alphabet)
            else:
                #Can we cope with a step as well? Problem is slice of slice...
                raise NotImplementedError
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
            s = Seq(self._load_seq(), self._alphabet)
            self._seq = s
            return s
    seq = property(fget=_get_seq,
                   fset=SeqRecord._set_seq, #Defined in base classes
                   doc="The sequence itself, as a Seq or MutableSeq object.")

    def _get_id(self):
        try:
            return self._id
        except AttributeError:
            #Load it now
            temp = self._load_id()
            self._id = temp
            return temp
    def _set_id(self, value):
        self._id = value
    id = property(fget=_get_id,
                  fset=_set_id,
                  doc="The sequence identifier (string)")

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
            temp = self._load_per_letter_annotations()
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
            temp = self._load_features()
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
            temp = self._load_dbxrefs()
            self._dbxrefs = temp
            return temp
    def _set_dbxrefs(self, value):
        if not isinstance(value, list):
            raise TypeError
        self._dbxrefs = value
    dbxrefs = property(fget=_get_features,
                       fset=_set_features,
                       doc="Database cross-references (list)")

    def _load_seq(self):
        """Extracts the (sub)sequence from the file."""
        raise NotImplementedError

    def _load_id(self):
        raise NotImplementedError

    def _load_name(self):
        return ""

    def _load_description(self):
        return ""

    def _load_annotations(self):
        return {}

    def _load_per_letter_annotations(self):
        return {}

    def _load_features(self):
        return []

    def _load_dbxrefs(self):
        return []
