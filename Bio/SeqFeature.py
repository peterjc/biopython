# Copyright 2000-2003 Jeff Chang.
# Copyright 2001-2008 Brad Chapman.
# Copyright 2005-2011 by Peter Cock.
# Copyright 2006-2009 Michiel de Hoon.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Represent a Sequence Feature holding info about a part of a sequence.

This is heavily modeled after the Biocorba SeqFeature objects, and
may be pretty biased towards GenBank stuff since I'm writing it
for the GenBank parser output...

What's here:

Base class to hold a Feature.
----------------------------
classes:
o SeqFeature

Hold information about a Reference.
----------------------------------

This is an attempt to create a General class to hold Reference type
information.

classes:
o Reference

Specify locations of a feature on a Sequence.
---------------------------------------------

This aims to handle, in Ewan's words, 'the dreaded fuzziness issue' in
much the same way as Biocorba. This has the advantages of allowing us
to handle fuzzy stuff in case anyone needs it, and also be compatible
with Biocorba.

classes:
o FeatureLocation - Specify the start and end location of a feature.
o CompoundLocation - Collection of FeatureLocation objects (for joins etc).

o ExactPosition - Specify the position as being exact.
o WithinPosition - Specify a position occuring within some range.
o BetweenPosition - Specify a position occuring between a range (OBSOLETE?).
o BeforePosition - Specify the position as being found before some base.
o AfterPosition - Specify the position as being found after some base.
o OneOfPosition - Specify a position where the location can be multiple positions.
"""

from Bio.Seq import MutableSeq, reverse_complement

class SeqFeature(object):
    """Represent a Sequence Feature on an object.

    Attributes:
    o location - the location of the feature on the sequence (FeatureLocation)
    o type - the specified type of the feature (ie. CDS, exon, repeat...)
    o location_operator - a string specifying how this SeqFeature may
    be related to others. For example, in the example GenBank feature
    shown below, the location_operator would be "join". This is a proxy
    for feature.location.operator and only applies to compound locations.
    o strand - A value specifying on which strand (of a DNA sequence, for
    instance) the feature deals with. 1 indicates the plus strand, -1 
    indicates the minus strand, 0 indicates stranded but unknown (? in GFF3),
    while the default of None indicates that strand doesn't apply (dot in GFF3,
    e.g. features on proteins). Note this is a shortcut for accessing the
    strand property of the feature's location.
    o id - A string identifier for the feature.
    o ref - A reference to another sequence. This could be an accession
    number for some different sequence. Note this is a shortcut for the
    reference property of the feature's location.
    o ref_db - A different database for the reference accession number.
    Note this is a shortcut for the reference property of the location
    o qualifiers - A dictionary of qualifiers on the feature. These are
    analagous to the qualifiers from a GenBank feature table. The keys of
    the dictionary are qualifier names, the values are the qualifier
    values.
    o sub_features - Obsolete attribute no longer used.

    When dealing with compound or split locations (e.g. joins in GenBank),
    Previous versions of Biopython would introduce additional SeqFeatures
    kept in the sub_features list of the 'parent' feature. For instance,
    if we having something like:

    CDS    join(1..10,30..40,50..60)

    Then the top level feature would be of type 'CDS' from 1 to 60 (actually 0
    to 60 in Python counting) with location_operator='join', and the three sub-
    features would also be of type 'CDS', and would be from 1 to 10, 30 to
    40 and 50 to 60, respectively (although actually using Python counting).

    These special locations are now handled with a special location object,
    a new CompoundLocation which holds a list of FeatureLocation objects.

    To get the nucleotide sequence for this CDS, you would need to take the
    parent sequence and do seq[0:10]+seq[29:40]+seq[49:60] (Python counting).
    Things are more complicated with strands and fuzzy positions. To save you
    dealing with all these special cases, the SeqFeature provides an extract
    method to do this for you.
    """
    def __init__(self, location = None, type = '', location_operator = '',
                 strand = None, id = "<unknown id>", 
                 qualifiers = None, sub_features = None,
                 ref = None, ref_db = None):
        """Initialize a SeqFeature on a Sequence.

        location can either be a FeatureLocation (with strand argument also
        given if required), or None.

        e.g. With no strand, on the forward strand, and on the reverse strand:

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> f1 = SeqFeature(FeatureLocation(5, 10), type="domain")
        >>> f1.strand == f1.location.strand == None
        True
        >>> f2 = SeqFeature(FeatureLocation(7, 110, strand=1), type="CDS")
        >>> f2.strand == f2.location.strand == +1
        True
        >>> f3 = SeqFeature(FeatureLocation(9, 108, strand=-1), type="CDS")
        >>> f3.strand == f3.location.strand == -1
        True

        An invalid strand will trigger an exception:

        >>> f4 = SeqFeature(FeatureLocation(50, 60), strand=2)
        Traceback (most recent call last):
           ...
        ValueError: Strand should be +1, -1, 0 or None, not 2

        Similarly if set via the FeatureLocation directly:

        >>> loc4 = FeatureLocation(50, 60, strand=2)
        Traceback (most recent call last):
           ...
        ValueError: Strand should be +1, -1, 0 or None, not 2

        For exact start/end positions, an integer can be used (as shown above)
        as shorthand for the ExactPosition object. For non-exact locations, the
        FeatureLocation must be specified via the appropriate position objects.

        Note that the strand, ref and ref_db arguments to the SeqFeature are
        now obsolete and will be deprecated in a future release (which will
        give warning messages) and later removed. Set them via the location
        object instead.

        Note that location_operator and sub_features arguments can no longer
        be used, instead do this via the CompoundLocation object.
        """
        if location is not None and not isinstance(location, FeatureLocation) \
        and not isinstance(location, CompoundLocation):
            raise TypeError("FeatureLocation, CompoundLocation (or None) required for the location")
        self.location = location

        self.type = type
        if location_operator:
            raise ValueError("Sorry, set location operator via the CompoundLocation now")
        if strand is not None:
            #TODO - Deprecation warning
            self.strand = strand
        self.id = id
        if qualifiers is None:
            qualifiers = {}
        self.qualifiers = qualifiers
        if sub_features is not None:
            raise ValueError("The sub_features argument and associated attribute "
                             "are no longer supported, use a CompoundLocation instead")
        if ref is not None:
            #TODO - Deprecation warning
            self.ref = ref
        if ref_db is not None:
            #TODO - Deprecation warning
            self.ref_db = ref_db

    def _get_sub_features(self):
        raise AttributeError("This will be a deprecation warning and return []")
        import warnings
        from Bio import BiopythonDeprecationWarning
        warnings.warn("The sub_features attribute is no longer used. Instead "
                      "see the SeqFeature's location is a CompoundLocation.",
                      BiopythonDeprecationWarning)
        return []
    sub_features = property(fget = _get_sub_features,
                            doc = "No longer used, see CompoundLocation object instead.")

    def _get_strand(self):
        return self.location.strand
    def _set_strand(self, value):
        try:
            self.location.strand = value
        except AttributeError:
            if self.location is None:
                if value is not None:
                    raise ValueError("Can't set strand without a location.")
            else:
                raise
    strand = property(fget = _get_strand, fset = _set_strand,
                      doc = """Feature's strand

                            This is a shortcut for feature.location.strand
                            """)

    def _get_ref(self):
        try:
            return self.location.ref
        except AttributeError:
            return None
    def _set_ref(self, value):
        try:
            self.location.ref = value
        except AttributeError:
            if self.location is None:
                if value is not None:
                    raise ValueError("Can't set ref without a location.")
            else:
                raise
    ref = property(fget = _get_ref, fset = _set_ref,
                   doc = """Feature location reference (e.g. accession).

                         This is a shortcut for feature.location.ref
                         """)

    def _get_ref_db(self):
        try:
            return self.location.ref_db
        except AttributeError:
            return None
    def _set_ref_db(self, value):
        self.location.ref_db = value
    ref_db = property(fget = _get_ref_db, fset = _set_ref_db,
                      doc = """Feature location reference's database.

                            This is a shortcut for feature.location.ref_db
                            """)

    def _get_location_operator(self):
        try:
            return self.location.operator
        except AttributeError:
            return None
    location_operator = property(fget = _get_location_operator,
                                 doc = "Location operator for compound locations (e.g. join).")

    def __repr__(self):
        """A string representation of the record for debugging."""
        answer = "%s(%s" % (self.__class__.__name__, repr(self.location))
        if self.type:
            answer += ", type=%s" % repr(self.type)
        if self.id and self.id != "<unknown id>":
            answer += ", id=%s" % repr(self.id)
        answer += ")"
        return answer

    def __str__(self):
        """A readable summary of the feature intended to be printed to screen.
        """
        out = "type: %s\n" % self.type
        out += "location: %s\n" % self.location
        if self.id and self.id != "<unknown id>":
            out += "id: %s\n" % self.id
        if self.ref or self.ref_db:
            out += "ref: %s:%s\n" % (self.ref, self.ref_db)
        out += "qualifiers: \n"
        for qual_key in sorted(self.qualifiers):
            out += "    Key: %s, Value: %s\n" % (qual_key,
                                               self.qualifiers[qual_key])
        return out

    def _shift(self, offset):
        """Returns a copy of the feature with its location shifted (PRIVATE).

        The annotation qaulifiers are copied."""
        return SeqFeature(location = self.location._shift(offset),
            type = self.type,
            id = self.id,
            qualifiers = dict(self.qualifiers.iteritems()))

    def _flip(self, length):
        """Returns a copy of the feature with its location flipped (PRIVATE).
        
        The argument length gives the length of the parent sequence. For
        example a location 0..20 (+1 strand) with parent length 30 becomes
        after flipping 10..30 (-1 strand). Strandless (None) or unknown
        strand (0) remain like that - just their end points are changed.

        The annotation qaulifiers are copied.
        """
        return SeqFeature(location = self.location._flip(length),
            type = self.type,
            id = self.id,
            qualifiers = dict(self.qualifiers.iteritems()))
    
    def extract(self, parent_sequence):
        """Extract feature sequence from the supplied parent sequence.

        The parent_sequence can be a Seq like object or a string, and will
        generally return an object of the same type. The exception to this is
        a MutableSeq as the parent sequence will return a Seq object.

        This should cope with complex locations including complements, joins
        and fuzzy positions. Even mixed strand features should work! This
        also covers features on protein sequences (e.g. domains), although
        here reverse strand features are not permitted.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_protein
        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> seq = Seq("MKQHKAMIVALIVICITAVVAAL", generic_protein)
        >>> f = SeqFeature(FeatureLocation(8,15), type="domain")
        >>> f.extract(seq)
        Seq('VALIVIC', ProteinAlphabet())

        Note - currently only compound locations of type "join" are supported.
        """
        return self.location.extract(parent_sequence)
    
    def __nonzero__(self):
        """Returns True regardless of the length of the feature.

        This behaviour is for backwards compatibility, since until the
        __len__ method was added, a SeqFeature always evaluated as True.

        Note that in comparison, Seq objects, strings, lists, etc, will all
        evaluate to False if they have length zero.

        WARNING: The SeqFeature may in future evaluate to False when its
        length is zero (in order to better match normal python behaviour)!
        """
        return True

    def __len__(self):
        """Returns the length of the region described by a feature.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_protein
        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> seq = Seq("MKQHKAMIVALIVICITAVVAAL", generic_protein)
        >>> f = SeqFeature(FeatureLocation(8,15), type="domain")
        >>> len(f)
        7
        >>> f.extract(seq)
        Seq('VALIVIC', ProteinAlphabet())
        >>> len(f.extract(seq))
        7

        For simple features without subfeatures this is the same as the region
        spanned (end position minus start position). However, for a feature
        defined by combining several subfeatures (e.g. a CDS as the join of
        several exons) the gaps are not counted (e.g. introns). This ensures
        that len(f) == len(f.extract(parent_seq)), and also makes sure things
        work properly with features wrapping the origin etc.
        """
        return len(self.location)

    def __iter__(self):
        """Iterate over the parent positions within the feature.

        The iteration order is strand aware, and can be thought of as moving
        along the feature using the parent sequence coordinates:

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> f = SeqFeature(FeatureLocation(5,10), type="domain", strand=-1)
        >>> len(f)
        5
        >>> for i in f: print i
        9
        8
        7
        6
        5
        >>> list(f)
        [9, 8, 7, 6, 5]
        """
        return iter(self.location)

    def __contains__(self, value):
        """Check if an integer position is within the feature.

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> f = SeqFeature(FeatureLocation(5,10), type="domain", strand=-1)
        >>> len(f)
        5
        >>> [i for i in range(15) if i in f]
        [5, 6, 7, 8, 9]

        For example, to see which features include a SNP position, you could
        use this:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("GenBank/NC_000932.gb", "gb")
        >>> for f in record.features:
        ...     if 1750 in f:
        ...         print f.type, f.strand, f.location
        source 1 [0:154478]
        gene -1 c[1716:4347]
        tRNA -1 join of c[1716:1751], c[4310:4347]

        Note that for a feature defined as a join of several subfeatures (e.g.
        the union of several exons) the gaps are not checked (e.g. introns).
        In this example, the tRNA location is defined in the GenBank file as
        complement(join(1717..1751,4311..4347)), so that position 1760 falls
        in the gap:

        >>> for f in record.features:
        ...     if 1760 in f:
        ...         print f.type, f.strand, f.location
        source 1 [0:154478]
        gene -1 c[1716:4347]

        Note that additional care may be required with fuzzy locations, for
        example just before a BeforePosition:

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> from Bio.SeqFeature import BeforePosition
        >>> f = SeqFeature(FeatureLocation(BeforePosition(3),8), type="domain")
        >>> len(f)
        5
        >>> [i for i in range(10) if i in f]
        [3, 4, 5, 6, 7]
        """
        if not isinstance(value, int):
            raise ValueError("Currently we only support checking for integer "
                             "positions being within a SeqFeature.")
        return value in self.location


# --- References

# TODO -- Will this hold PubMed and Medline information decently?
class Reference(object):
    """Represent a Generic Reference object.

    Attributes:
    o location - A list of Location objects specifying regions of
    the sequence that the references correspond to. If no locations are
    specified, the entire sequence is assumed.
    o authors - A big old string, or a list split by author, of authors
    for the reference.
    o title - The title of the reference.
    o journal - Journal the reference was published in.
    o medline_id - A medline reference for the article.
    o pubmed_id - A pubmed reference for the article.
    o comment - A place to stick any comments about the reference.
    """
    def __init__(self):
        self.location = []
        self.authors = ''
        self.consrtm = ''
        self.title = ''
        self.journal = ''
        self.medline_id = ''
        self.pubmed_id = ''
        self.comment = ''

    def __str__(self):
        """Output an informative string for debugging.
        """
        out = ""
        for single_location in self.location:
            out += "location: %s\n" % single_location
        out += "authors: %s\n" % self.authors
        if self.consrtm:
            out += "consrtm: %s\n" % self.consrtm
        out += "title: %s\n" % self.title
        out += "journal: %s\n" % self.journal
        out += "medline id: %s\n" % self.medline_id
        out += "pubmed id: %s\n" % self.pubmed_id
        out += "comment: %s\n" % self.comment
        return out

    def __repr__(self):
        #TODO - Update this is __init__ later accpets values
        return "%s(title=%s, ...)" % (self.__class__.__name__,
                                      repr(self.title))

# --- Handling feature locations

class FeatureLocation(object):
    """Specify the location of a feature along a sequence.

    This attempts to deal with fuzziness of position ends, but also
    make it easy to get the start and end in the 'normal' case (no
    fuzziness).

    You should access the start and end attributes with
    your_location.start and your_location.end. If the start and
    end are exact, this will return the positions, if not, we'll return
    the approriate Fuzzy class with info about the position and fuzziness.

    Note that the start and end location numbering follow Python's scheme,
    thus a GenBank entry of 123..150 (one based counting) becomes a location
    of [122:150] (zero based counting).
    """
    def __init__(self, start, end, strand=None, ref=None, ref_db=None):
        """Specify the start, end, strand etc of a sequence feature.

        start and end arguments specify the values where the feature begins
        and ends. These can either by any of the *Position objects that
        inherit from AbstractPosition, or can just be integers specifying the
        position. In the case of integers, the values are assumed to be
        exact and are converted in ExactPosition arguments. This is meant
        to make it easy to deal with non-fuzzy ends.

        i.e. Short form:
        
        >>> from Bio.SeqFeature import FeatureLocation
        >>> loc = FeatureLocation(5,10)
        
        Explicit form:

        >>> from Bio.SeqFeature import FeatureLocation, ExactPosition
        >>> loc = FeatureLocation(ExactPosition(5),ExactPosition(10))

        Other fuzzy positions are used similarly,

        >>> from Bio.SeqFeature import FeatureLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc2 = FeatureLocation(BeforePosition(5),AfterPosition(10))

        For nucleotide features you will also want to specify the strand,
        use 1 for the forward (plus) strand, -1 for the reverse (negative)
        strand, 0 for stranded but strand unknown (? in GFF3), or None for
        when the strand does not apply (dot in GFF3), e.g. features on
        proteins.

        >>> loc = FeatureLocation(5, 10, strand=+1)
        >>> print loc.strand
        1

        Normally feature locations are given relative to the parent
        sequence you are working with, but an explicit accession can
        be given with the optional ref and db_ref strings:

        >>> loc = FeatureLocation(105172, 108462, ref="AL391218.9")
        >>> print loc.ref
        AL391218.9

        """
        if isinstance(start, AbstractPosition):
            self._start = start
        else:
            self._start = ExactPosition(start)
        if isinstance(end, AbstractPosition):
            self._end = end
        else:
            self._end = ExactPosition(end)
        self.strand = strand
        self.ref = ref
        self.ref_db = ref_db

    def _get_strand(self):
        return self._strand
    def _set_strand(self, value):
        if value not in [+1, -1, 0, None]:
            raise ValueError("Strand should be +1, -1, 0 or None, not %r" \
                             % value)
        self._strand = value
    strand = property(fget = _get_strand, fset = _set_strand,
                      doc = "Strand of the location (+1, -1, 0 or None).")

    def __str__(self):
        """Returns a representation of the location (with python counting).

        For the simple case this uses the python splicing syntax, [122:150]
        (zero based counting) which GenBank would call 123..150 (one based
        counting).
        """
        if self.strand == -1:
            return "c[%s:%s]" % (self._start, self._end)
        return "[%s:%s]" % (self._start, self._end)

    def __repr__(self):
        """A string representation of the location for debugging."""
        optional = ""
        if self.strand is not None:
            optional += ", strand=%r" % self.strand
        if self.ref is not None:
            optional += ", ref=%r" % self.ref
        if self.ref_db is not None:
            optional += ", ref_db=%r" % self.ref_db
        return "%s(%r, %r%s)" \
                   % (self.__class__.__name__, self.start, self.end, optional)

    def __nonzero__(self):
        """Returns True regardless of the length of the feature.

        This behaviour is for backwards compatibility, since until the
        __len__ method was added, a FeatureLocation always evaluated as True.

        Note that in comparison, Seq objects, strings, lists, etc, will all
        evaluate to False if they have length zero.

        WARNING: The FeatureLocation may in future evaluate to False when its
        length is zero (in order to better match normal python behaviour)!
        """
        return True

    def __len__(self):
        """Returns the length of the region described by the FeatureLocation.
        
        Note that extra care may be needed for fuzzy locations, e.g.

        >>> from Bio.SeqFeature import FeatureLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc = FeatureLocation(BeforePosition(5),AfterPosition(10))
        >>> len(loc)
        5
        """
        #TODO - Should we use nofuzzy_start and nofuzzy_end here?
        return self._end.position + self._end.extension - self._start.position

    def __contains__(self, value):
        """Check if an integer position is within the FeatureLocation.

        Note that extra care may be needed for fuzzy locations, e.g.

        >>> from Bio.SeqFeature import FeatureLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc = FeatureLocation(BeforePosition(5),AfterPosition(10))
        >>> len(loc)
        5
        >>> [i for i in range(15) if i in loc]
        [5, 6, 7, 8, 9]
        """
        if not isinstance(value, int):
            raise ValueError("Currently we only support checking for integer "
                             "positions being within a FeatureLocation.")
        #TODO - Should we use nofuzzy_start and nofuzzy_end here?
        if value < self._start.position \
        or value >= self._end.position + self._end.extension:
            return False
        else:
            return True

    def __iter__(self):
        """Iterate over the parent positions within the FeatureLocation.

        >>> from Bio.SeqFeature import FeatureLocation
        >>> from Bio.SeqFeature import BeforePosition, AfterPosition
        >>> loc = FeatureLocation(BeforePosition(5),AfterPosition(10))
        >>> len(loc)
        5
        >>> for i in loc: print i
        5
        6
        7
        8
        9
        >>> list(loc)
        [5, 6, 7, 8, 9]
        >>> [i for i in range(15) if i in loc]
        [5, 6, 7, 8, 9]

        Note this is strand aware:

        >>> loc = FeatureLocation(BeforePosition(5), AfterPosition(10), strand = -1)
        >>> list(loc)
        [9, 8, 7, 6, 5]
        """
        if self.strand == -1:
            return iter(xrange(self.nofuzzy_end-1, self.nofuzzy_start-1, -1))
        else:
            return iter(xrange(self.nofuzzy_start, self.nofuzzy_end))

    def _shift(self, offset):
        """Returns a copy of the location shifted by the offset (PRIVATE)."""
        return FeatureLocation(start = self._start._shift(offset),
                               end = self._end._shift(offset),
                               strand = self.strand)

    def _flip(self, length):
        """Returns a copy of the location after the parent is reversed (PRIVATE)."""
        #Note this will flip the start and end too!
        if self.strand == +1:
            flip_strand = -1
        elif self.strand == -1:
            flip_strand = +1
        else:
            #0 or None
            flip_strand = self.strand
        return FeatureLocation(start = self._end._flip(length),
                               end = self._start._flip(length),
                               strand = flip_strand)

    start = property(fget= lambda self : self._start,
                 doc="Start location (possibly a fuzzy position, read only).")

    end = property(fget= lambda self : self._end,
                   doc="End location (possibly a fuzzy position, read only).")

    nofuzzy_start = property(
        fget=lambda self: self._start.position,
        doc="""Start position (integer, approximated if fuzzy, read only).

        To get non-fuzzy attributes (ie. the position only) ask for
        'location.nofuzzy_start', 'location.nofuzzy_end'. These should return
        the largest range of the fuzzy position. So something like:
        (10.20)..(30.40) should return 10 for start, and 40 for end.
        """)

    nofuzzy_end = property(
        fget=lambda self: self._end.position + self._end.extension,
        doc="""End position (integer, approximated if fuzzy, read only).

        To get non-fuzzy attributes (ie. the position only) ask for
        'location.nofuzzy_start', 'location.nofuzzy_end'. These should return
        the largest range of the fuzzy position. So something like:
        (10.20)..(30.40) should return 10 for start, and 40 for end.
        """)

    def extract(self, parent_sequence):
        """Extract feature sequence from the supplied parent sequence."""
        if isinstance(parent_sequence, MutableSeq):
            #This avoids complications with reverse complements
            #(the MutableSeq reverse complement acts in situ)
            parent_sequence = parent_sequence.toseq()
        f_seq = parent_sequence[self.nofuzzy_start:\
                                self.nofuzzy_end]
        if self.strand == -1:
            #TODO - MutableSeq?
            try:
                f_seq = f_seq.reverse_complement()
            except AttributeError:
                assert isinstance(f_seq, str)
                f_seq = reverse_complement(f_seq)
        return f_seq


class CompoundLocation(object):
    """For handling joins etc where a feature location has several parts."""
    def __init__(self, parts, operator="join"):
        """Create a compound location with several parts.

        >>> from Bio.SeqFeature import FeatureLocation, CompoundLocation
        >>> f1 = FeatureLocation(10, 40, strand=+1)
        >>> f2 = FeatureLocation(50, 59, strand=+1)
        >>> f = CompoundLocation([f1, f2])
        >>> len(f) == len(f1) + len(f2) == 39 == len(list(f))
        True
        >>> print f.operator
        join
        >>> 5 in f
        False
        >>> 15 in f
        True
        >>> f.strand
        1

        Notice that the strand of the compound location is computed
        automatically - in the case of mixed strands on the sub-locations
        the overall strand is set to None.

        >>> f = CompoundLocation([FeatureLocation(3, 6, strand=+1),
        ...                       FeatureLocation(10, 13, strand=-1)])
        >>> print f.strand
        None
        >>> len(f)
        6
        >>> list(f)
        [3, 4, 5, 12, 11, 10]
        """
        self.operator = operator
        self.parts = list(parts)
        for loc in self.parts:
            if not isinstance(loc, FeatureLocation):
                raise ValueError("CompoundLocation should be given a list of "
                                 "FeatureLocation objects, not %s" % loc.__class__)
        if len(self.parts) < 2:
            raise ValueError("CompoundLocation should have at least 2 parts - have %i" % len(self.parts))

    def __str__(self):
        """Returns a representation of the location (with python counting)."""
        return "%s of %s" % (self.operator, ", ".join(str(loc) for loc in self.parts))

    def __repr__(self):
        """String representation of the location for debugging."""
        return "%s(%r, %r)" % (self.__class__.__name__, \
                               self.parts, self.operator)

    def _get_strand(self):
        # Historically a join on the reverse strand has been represented
        # in Biopython with both the parent SeqFeature and its children
        # (the exons for a CDS) all given a strand of -1.  Likewise, for
        # a join feature on the forward strand they all have strand +1.
        # However, we must also consider evil mixed strand examples like
        # this, join(complement(69611..69724),139856..140087,140625..140650)
        if len(set(loc.strand for loc in self.parts))==1:
            return self.parts[0].strand
        else:
            return None # i.e. mixed strands
    strand = property(fget = _get_strand,
                      doc = "Overall strand of the compound location (read only).")

    def __contains__(self, value):
        """Check if an integer position is within the location."""
        for loc in self.parts:
            if value in loc:
                return True
        return False

    def __nonzero__(self):
        """Returns True regardless of the length of the feature.

        This behaviour is for backwards compatibility, since until the
        __len__ method was added, a FeatureLocation always evaluated as True.

        Note that in comparison, Seq objects, strings, lists, etc, will all
        evaluate to False if they have length zero.

        WARNING: The FeatureLocation may in future evaluate to False when its
        length is zero (in order to better match normal python behaviour)!
        """
        return True

    def __len__(self):
        return sum(len(loc) for loc in self.parts)

    def __iter__(self):
        #TODO - Fix GenBank parser to avoid this hack
        if self.strand == -1:
            for loc in self.parts:
                assert loc.strand == -1
            for loc in self.parts[::-1]:
                for pos in loc:
                    yield pos
        else:
            for loc in self.parts:
                for pos in loc:
                    yield pos

    def _shift(self, offset):
        """Returns a copy of the location shifted by the offset (PRIVATE)."""
        return CompoundLocation([loc._shift(offset) for loc in self.parts],
                                self.operator)

    def _flip(self, length):
        """Returns a copy of the location after the parent is reversed (PRIVATE).

        Note that the order of the parts is reversed too.
        """
        return CompoundLocation([loc._flip(length) for loc in self.parts[::-1]],
                                self.operator)

    def _get_start(self):
        starts = set(loc.start for loc in self.parts)
        smallest = min(pos.position for pos in starts)
        starts = [pos for pos in starts if pos.position == smallest]
        #Due to fuzzyness, could in theory have more than one
        #(overlapping) sub-part with same start...
        return starts[0]
    start = property(fget=_get_start,
        doc = """Start position (integer, approximated if fuzzy, read only).

        See comments for the start attribute.
        """)

    nofuzzy_start = property(
        fget=lambda self: min(loc.start.position for loc in self.parts),
        doc="""Smallest start position (object) of sub-parts (read only).

        Note that the start/end are intended to give the numerical range
        spanned by this feature suitable for drawing etc. In the case of
        a circular genome and a feature spanning the origin represented
        by two parts, this means the CompoundLocation's start will be 0
        and its end will be the length of the genome. 
        """)

    def _get_end(self):
        ends = set(loc.end for loc in self.parts)
        biggest = max(pos.position + pos.extension for pos in ends)
        #Due to fuzzyness, could in theory have more than one
        #(overlapping) sub-part with same end...
        ends = [pos for pos in ends if pos.position + pos.extension == biggest]
        return ends[0]
    end = property(
        fget=_get_end,
        doc="""Largest end position (object) of sub-parts (read only).

        Note that the start/end are intended to give the numerical range
        spanned by this feature suitable for drawing etc. In the case of
        a circular genome and a feature spanning the origin represented
        by two parts, this means the CompoundLocation's start will be 0
        and its end will be the length of the genome.
        """)

    nofuzzy_end = property(
        fget=lambda self: max(loc.end.position + loc.end.extension for loc in self.parts),
        doc="""End position (integer, approximated if fuzzy, read only).

        See the comments for the end position attribute.
        """)

    def extract(self, parent_sequence):
        """Extract feature sequence from the supplied parent sequence."""
        if isinstance(parent_sequence, MutableSeq):
            #This avoids complications with reverse complements
            #(the MutableSeq reverse complement acts in situ)
            parent_sequence = parent_sequence.toseq()
        if self.operator!="join":
            raise ValueError(self.operator)
        if self.strand == -1:
            #This is a special case given how the GenBank parser works.
            #Must avoid doing the reverse complement twice.
            #TODO - Fix this as part of sub_feature -> CompoundLocation
            parts = []
            for loc in self.parts:
                assert loc.strand==-1
                parts.append(parent_sequence[loc.nofuzzy_start:\
                                             loc.nofuzzy_end])
            f_seq = parts[0]
            for part in parts[1:]:
                f_seq += part
            try:
                f_seq = f_seq.reverse_complement()
            except AttributeError:
                assert isinstance(f_seq, str)
                f_seq = reverse_complement(f_seq)
        else:
            #This copes with mixed strand features:
            parts = [loc.extract(parent_sequence) for loc in self.parts]
            #We use addition rather than a join to avoid alphabet issues:
            f_seq = parts[0]
            for part in parts[1:]:
                f_seq += part
        return f_seq


class AbstractPosition(object):
    """Abstract base class representing a position.
    """
    def __init__(self, position, extension):
        self.position = position
        assert extension >= 0, extension
        self.extension = extension

    def __repr__(self):
        """String representation of the location for debugging."""
        return "%s(%s,%s)" % (self.__class__.__name__, \
                              repr(self.position), repr(self.extension))

    def __hash__(self):
        """Simple position based hash."""
        #Note __hash__ must be implemented on Python 3.x if overriding __eq__
        return hash(self.position)

    def __eq__(self, other):
        """A simple equality for positions.

        This is very simple-minded and just compares the position attribute
        of the features; extensions are not considered at all. This could
        potentially be expanded to try to take advantage of extensions.
        """
        assert isinstance(other, AbstractPosition), \
          "We can only do comparisons between Biopython Position objects."
        return self.position == other.position

    def __ne__(self, other):
        """A simple non-equality for positions.

        This is very simple-minded and just compares the position attribute
        of the features; extensions are not considered at all. This could
        potentially be expanded to try to take advantage of extensions.
        """
        assert isinstance(other, AbstractPosition), \
          "We can only do comparisons between Biopython Position objects."
        return self.position != other.position

    def __le__(self, other):
        """A simple less than or equal for positions.

        This is very simple-minded and just compares the position attribute
        of the features; extensions are not considered at all. This could
        potentially be expanded to try to take advantage of extensions.
        """
        assert isinstance(other, AbstractPosition), \
          "We can only do comparisons between Biopython Position objects."
        return self.position <= other.position

    def __lt__(self, other):
        """A simple less than or equal for positions.

        This is very simple-minded and just compares the position attribute
        of the features; extensions are not considered at all. This could
        potentially be expanded to try to take advantage of extensions.
        """
        assert isinstance(other, AbstractPosition), \
          "We can only do comparisons between Biopython Position objects."
        return self.position < other.position

    def __ge__(self, other):
        """A simple less than or equal for positions.

        This is very simple-minded and just compares the position attribute
        of the features; extensions are not considered at all. This could
        potentially be expanded to try to take advantage of extensions.
        """
        assert isinstance(other, AbstractPosition), \
          "We can only do comparisons between Biopython Position objects."
        return self.position >= other.position

    def __gt__(self, other):
        """A simple less than or equal for positions.

        This is very simple-minded and just compares the position attribute
        of the features; extensions are not considered at all. This could
        potentially be expanded to try to take advantage of extensions.
        """
        assert isinstance(other, AbstractPosition), \
          "We can only do comparisons between Biopython Position objects."
        return self.position > other.position

    def _shift(self, offset):
        #We want this to maintain the subclass when called from a subclass
        return self.__class__(self.position + offset, self.extension)

    def _flip(self, length):
        #We want this to maintain the subclass when called from a subclass
        return self.__class__(length - self.position - self.extension,
                              self.extension)


class ExactPosition(AbstractPosition):
    """Specify the specific position of a boundary.

    o position - The position of the boundary.
    o extension - An optional argument which must be zero since we don't
    have an extension. The argument is provided so that the same number of
    arguments can be passed to all position types.

    In this case, there is no fuzziness associated with the position.
    """
    def __init__(self, position, extension = 0):
        if extension != 0:
            raise AttributeError("Non-zero extension %s for exact position."
                                 % extension)
        AbstractPosition.__init__(self, position, 0)

    def __repr__(self):
        """String representation of the ExactPosition location for debugging."""
        assert self.extension == 0
        return "%s(%s)" % (self.__class__.__name__, repr(self.position))

    def __str__(self):
        return str(self.position)

class UncertainPosition(ExactPosition):
    """Specify a specific position which is uncertain.
    
    This is used in UniProt, e.g. ?222 for uncertain position 222, or in the
    XML format explicitly marked as uncertain. Does not apply to GenBank/EMBL.
    """
    pass

class UnknownPosition(AbstractPosition):
    """Specify a specific position which is unknown (has no position).

    This is used in UniProt, e.g. ? or in the XML as unknown.
    """
    def __init__(self):
        self.position = None
        self.extension = None
        pass

    def __repr__(self):
        """String representation of the UnknownPosition location for debugging."""
        return "%s()" % self.__class__.__name__
        
class WithinPosition(AbstractPosition):
    """Specify the position of a boundary within some coordinates.

    Arguments:
    o position - The start position of the boundary
    o extension - The range to which the boundary can extend.

    This allows dealing with a position like ((1.4)..100). This
    indicates that the start of the sequence is somewhere between 1
    and 4. To represent that with this class we would set position as
    1 and extension as 3.
    """
    def __init__(self, position, extension = 0):
        AbstractPosition.__init__(self, position, extension)

    def __str__(self):
        return "(%s.%s)" % (self.position, self.position + self.extension)


class BetweenPosition(AbstractPosition):
    """Specify the position of a boundary between two coordinates (OBSOLETE?).

    Arguments:
    o position - The start position of the boundary.
    o extension - The range to the other position of a boundary.

    This specifies a coordinate which is found between the two positions.
    So this allows us to deal with a position like ((1^2)..100). To
    represent that with this class we set position as 1 and the
    extension as 1.
    """
    def __init__(self, position, extension = 0):
        AbstractPosition.__init__(self, position, extension)

    def __str__(self):
        return "(%s^%s)" % (self.position, self.position + self.extension)


class BeforePosition(AbstractPosition):
    """Specify a position where the actual location occurs before it.

    Arguments:
    o position - The upper boundary of where the location can occur.
    o extension - An optional argument which must be zero since we don't
    have an extension. The argument is provided so that the same number of
    arguments can be passed to all position types.

    This is used to specify positions like (<10..100) where the location
    occurs somewhere before position 10.
    """
    def __init__(self, position, extension = 0):
        if extension != 0:
            raise AttributeError("Non-zero extension %s for exact position."
                                 % extension)
        AbstractPosition.__init__(self, position, 0)

    def __repr__(self):
        """A string representation of the location for debugging."""
        assert self.extension == 0
        return "%s(%s)" % (self.__class__.__name__, repr(self.position))

    def __str__(self):
        return "<%s" % self.position

    def _flip(self, length):
        return AfterPosition(length - self.position)

class AfterPosition(AbstractPosition):
    """Specify a position where the actual location is found after it.

    Arguments:
    o position - The lower boundary of where the location can occur.
    o extension - An optional argument which must be zero since we don't
    have an extension. The argument is provided so that the same number of
    arguments can be passed to all position types.

    This is used to specify positions like (>10..100) where the location
    occurs somewhere after position 10.
    """
    def __init__(self, position, extension = 0):
        if extension != 0:
            raise AttributeError("Non-zero extension %s for exact position."
                                 % extension)
        AbstractPosition.__init__(self, position, 0)

    def __repr__(self):
        """A string representation of the location for debugging."""
        assert self.extension == 0
        return "%s(%s)" % (self.__class__.__name__, repr(self.position))

    def __str__(self):
        return ">%s" % self.position

    def _flip(self, length):
        return BeforePosition(length - self.position)


class OneOfPosition(AbstractPosition):
    """Specify a position where the location can be multiple positions.

    This models the GenBank 'one-of(1888,1901)' function, and tries
    to make this fit within the Biopython Position models. In our case
    the position of the "one-of" is set as the lowest choice, and the
    extension is the range to the highest choice.
    """
    def __init__(self, position_list):
        """Initialize with a set of posssible positions.

        position_list is a list of AbstractPosition derived objects,
        specifying possible locations.
        """
        # unique attribute for this type of positions
        self.position_choices = position_list
        # find the smallest and largest position in the choices
        smallest = None
        largest = None
        for position_choice in self.position_choices:
            assert isinstance(position_choice, AbstractPosition), \
              "Expected position objects, got %r" % position_choice
            if smallest is None and largest is None:
                smallest = position_choice.position
                largest = position_choice.position
            elif position_choice.position > largest:
                largest = position_choice.position
            elif position_choice.position < smallest:
                smallest = position_choice.position
        # initialize with our definition of position and extension
        AbstractPosition.__init__(self, smallest, largest - smallest)

    def __repr__(self):
        """String representation of the OneOfPosition location for debugging."""
        return "%s(%s)" % (self.__class__.__name__, \
                           repr(self.position_choices))

    def __str__(self):
        out = "one-of("
        for position in self.position_choices:
            out += "%s," % position
        # replace the last comma with the closing parenthesis
        out = out[:-1] + ")"
        return out

    def _shift(self, offset):
        return self.__class__([position_choice._shift(offset) \
                               for position_choice in self.position_choices])

    def _flip(self, length):
        return OneOfPosition([p._flip(length) for p in self.position_choices[::-1]])


class PositionGap(object):
    """Simple class to hold information about a gap between positions.
    """
    def __init__(self, gap_size):
        """Intialize with a position object containing the gap information.
        """
        self.gap_size = gap_size

    def __repr__(self):
        """A string representation of the position gap for debugging."""
        return "%s(%s)" % (self.__class__.__name__, repr(self.gap_size))
    
    def __str__(self):
        out = "gap(%s)" % self.gap_size
        return out

def _test():
    """Run the Bio.SeqFeature module's doctests (PRIVATE).

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..","Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
    elif os.path.isdir(os.path.join("Tests")) :
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("Tests"))
        doctest.testmod()
        os.chdir(cur_dir)
        del cur_dir
        print "Done"


if __name__ == "__main__":
    _test()
