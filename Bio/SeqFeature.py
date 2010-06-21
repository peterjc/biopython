# Copyright 2000-2003 Jeff Chang.
# Copyright 2001-2008 Brad Chapman.
# Copyright 2005-2010 by Peter Cock.
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
    shown below, the location_operator would be "join"
    o strand - A value specifying on which strand (of a DNA sequence, for
    instance) the feature deals with. 1 indicates the plus strand, -1 
    indicates the minus strand, 0 indicates both strands, and None indicates
    that strand doesn't apply (ie. for proteins) or is not known.
    o id - A string identifier for the feature.
    o ref - A reference to another sequence. This could be an accession
    number for some different sequence.
    o ref_db - A different database for the reference accession number.
    o qualifiers - A dictionary of qualifiers on the feature. These are
    analagous to the qualifiers from a GenBank feature table. The keys of
    the dictionary are qualifier names, the values are the qualifier
    values.
    o sub_features - Additional SeqFeatures which fall under this 'parent'
    feature. For instance, if we having something like:

    CDS    join(1..10,30..40,50..60)

    Then the top level feature would be of type 'CDS' from 1 to 60 (actually 0
    to 60 in Python counting) with location_operator='join', and the three sub-
    features would also be of type 'CDS', and would be from 1 to 10, 30 to
    40 and 50 to 60, respectively (although actually using Python counting).

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
        given if required), or a Python slice (with strand given as the step).

        e.g. With no strand, on the forward strand, and on the reverse strand:

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> f1 = SeqFeature(FeatureLocation(5,10), type="domain")
        >>> f2 = SeqFeature(FeatureLocation(7,110), strand=1, type="CDS")
        >>> f3 = SeqFeature(FeatureLocation(9,108), strand=-1, type="CDS")

        An invalid strand will trigger an exception:

        >>> f4 = SeqFeature(FeatureLocation(50,60), strand=2)
        Traceback (most recent call last):
           ...
        ValueError: Strand should be +1, -1, 0 or None, not 2

        For exact start/end positions, an integer can be used (as shown above)
        as shorthand for the ExactPosition object. For non-exact locations, the
        FeatureLocation must be specified via the appropriate position objects.
        """
        if strand not in [-1, 0, 1, None] :
            raise ValueError("Strand should be +1, -1, 0 or None, not %s" \
                             % repr(strand))
        if location is not None and not isinstance(location, FeatureLocation):
            raise TypeError("FeatureLocation (or None) required for the location")
        self.location = location

        self.type = type
        self.location_operator = location_operator
        self.strand = strand
        self.id = id
        if qualifiers is None:
            qualifiers = {}
        self.qualifiers = qualifiers
        if sub_features is None:
            sub_features = []
        self.sub_features = sub_features
        self.ref = ref 
        self.ref_db = ref_db

    def __repr__(self):
        """A string representation of the record for debugging."""
        answer = "%s(%s" % (self.__class__.__name__, repr(self.location))
        if self.type:
            answer += ", type=%s" % repr(self.type)
        if self.location_operator:
            answer += ", location_operator=%s" % repr(self.location_operator)
        if self.strand:
            answer += ", strand=%s" % repr(self.strand)
        if self.id and self.id != "<unknown id>":
            answer += ", id=%s" % repr(self.id)
        if self.ref:
            answer += ", ref=%s" % repr(self.ref)
        if self.ref_db:
            answer += ", ref_db=%s" % repr(self.ref_db)
        answer += ")"
        return answer

    def __str__(self):
        """A readable summary of the feature intended to be printed to screen.
        """
        out = "type: %s\n" % self.type
        out += "location: %s\n" % self.location
        if self.id:
            out += "id: %s\n" % self.id
        if self.ref or self.ref_db:
            out += "ref: %s:%s\n" % (self.ref, self.ref_db)
        out += "strand: %s\n" % self.strand
        out += "qualifiers: \n"
        qualifier_keys = self.qualifiers.keys()
        qualifier_keys.sort()
        for qual_key in qualifier_keys:
            out += "    Key: %s, Value: %s\n" % (qual_key,
                                               self.qualifiers[qual_key])
        if len(self.sub_features) != 0:
            out += "Sub-Features\n"
            for sub_feature in self.sub_features:
                out +="%s\n" % sub_feature

        return out

    def _shift(self, offset):
        """Returns a copy of the feature with its location shifted (PRIVATE).

        The annotation qaulifiers are copied."""
        return SeqFeature(location = self.location._shift(offset),
                type = self.type,
                location_operator = self.location_operator,
                strand = self.strand,
                id = self.id,
                qualifiers = dict(self.qualifiers.iteritems()),
                sub_features = [f._shift(offset) for f in self.sub_features],
                ref = self.ref,
                ref_db = self.ref_db)

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

        Note - currently only sub-features of type "join" are supported.
        """
        if isinstance(parent_sequence, MutableSeq):
            #This avoids complications with reverse complements
            #(the MutableSeq reverse complement acts in situ)
            parent_sequence = parent_sequence.toseq()
        if self.sub_features:
            if self.location_operator!="join":
                raise ValueError(self.location_operator)
            if self.strand == -1:
                #This is a special case given how the GenBank parser works.
                #Must avoid doing the reverse complement twice.
                parts = []
                for f_sub in self.sub_features:
                    assert f_sub.strand==-1
                    parts.append(parent_sequence[f_sub.location.nofuzzy_start:\
                                                 f_sub.location.nofuzzy_end])
            else:
                #This copes with mixed strand features:
                parts = [f_sub.extract(parent_sequence) \
                         for f_sub in self.sub_features]
            #We use addition rather than a join to avoid alphabet issues:
            f_seq = parts[0]
            for part in parts[1:] : f_seq += part
        else:
            f_seq = parent_sequence[self.location.nofuzzy_start:\
                                    self.location.nofuzzy_end]
        if self.strand == -1:
            #TODO - MutableSeq?
            try:
                f_seq = f_seq.reverse_complement()
            except AttributeError:
                assert isinstance(f_seq, str)
                f_seq = reverse_complement(f_seq)
        return f_seq
    
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
        if self.sub_features:
            return sum(len(f) for f in self.sub_features)
        else:
            return len(self.location)
    
    def get_local_coord(self, parent_coordinate):
        """Give the local coordinate for a parent coordinate.
        
        Input value is an integer (from 0 to the length of the parent), which
        will be converted to the features's coordindate system (giving an
        integer from 0 to the length of the feature). For the reverse mapping
        use the get_parent_coord method.
        
        See also the "in" functionality, a local coordinate can only be found
        if the position is in the feature:

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> f = SeqFeature(FeatureLocation(8,98), type="CDS", strand=+1)
        >>> len(f)
        90
        >>> 10 in f
        True
        >>> f.get_local_coord(10)
        2
        >>> f.get_parent_coord(2)
        10
        >>> 5 in f
        False
        >>> f.get_local_coord(5)
        Traceback (most recent call last):
            ...
        IndexError: Position 5 not in feature

        To demonstrate the two coordindate systems:

        >>> for local_coord in range(len(f)):
        ...    parent_coord = f.get_parent_coord(local_coord)
        ...    assert local_coord == f.get_local_coord(parent_coord)

        Here is an example on the reverse stand:

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> f = SeqFeature(FeatureLocation(10,100), type="CDS", strand=-1)
        >>> len(f)
        90
        >>> 99 in f
        True
        >>> f.get_local_coord(99)
        0
        >>> f.get_parent_coord(0)
        99
        >>> 100 in f
        False
        >>> 100 in f.location
        False
        >>> f.get_local_coord(100)
        Traceback (most recent call last):
            ...
        IndexError: Position 100 not in feature

        And again, to demonstrate the two coordindate systems:

        >>> for local_coord in range(len(f)):
        ...    parent_coord = f.get_parent_coord(local_coord)
        ...    assert local_coord == f.get_local_coord(parent_coord)

        """
        if parent_coordinate < 0:
            #TODO - Can we handle this case nicely to mean counting from the
            #end of the parent sequence? We'd need to know the length of the
            #parent sequence.
            raise ValueError("Parent co-ordinate should be at least zero")
        if self.sub_features:
            if self.strand == -1:
                #Reverse strand
                for f in self.sub_features[::-1]:
                    try:
                        local = f.get_local_coordinate(parent_coordinate)
                    except IndexError:
                        pass
            else:
                for f in self.sub_features:
                    try:
                        local = f.get_local_coordinate(parent_coordinate)
                    except IndexError:
                        pass
            assert False, "Problem getting co-ordinate via sub-features"
        elif parent_coordinate not in self.location:
            raise IndexError("Position %i not in feature" % parent_coordinate)
        elif self.strand == -1:
            return self.location.end.position + self.location.end.extension \
                   - parent_coordinate - 1
        else:
            return parent_coordinate - self.location.start.position
        
    def get_parent_coord(self, local_coordinate):
        """Give the parental coordinate for a local coordinate.
        
        Input value is an integer (from 0 to the length of the feature), which
        will be converted to the parental sequence's coordindate system. For
        the reverse mapping use the get_local_coord method.
        
        This is not intended to help you slice the parent sequence (due to the
        complications discussed below) but rather for accessing individual
        letters from the parent sequence such as SNPs.

        For a simple feature on the forward strand, this just adds the start
        position of the feature:

        >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
        >>> f = SeqFeature(FeatureLocation(8,98), type="CDS", strand=+1)
        >>> len(f)
        90
        >>> f.get_parent_coord(0)
        8
        >>> f.get_parent_coord(89)
        97
        
        For a simple feature on the reverse strand, this value is calculated
        relative to the end position (which is the biological start):

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_protein
        >>> f = SeqFeature(FeatureLocation(8,98), type="CDS", strand=-1)
        >>> len(f)
        90
        >>> f.get_parent_coord(0)
        97
        >>> f.get_parent_coord(1)
        96
        
        Let's look at a complicated CDS feature with a GenBank location string
        join(3462..3615,3698..3978,4077..4307,4408..4797,4876..5028,5141..5332)

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("GenBank/arab1.gb", "gb")
        >>> f = record.features[1]
        >>> print f.type, f.strand, len(f.sub_features)
        CDS 1 6
        >>> for sf in f.sub_features:
        ...    print sf.location, len(sf)
        [3461:3615] 154
        [3697:3978] 281
        [4076:4307] 231
        [4407:4797] 390
        [4875:5028] 153
        [5140:5332] 192
        >>> print f.location, len(f)
        [3461:5332] 1401
        
        Recall you can get the sequence of the feature with the extract method:

        >>> feature_seq = f.extract(record.seq)
        >>> len(feature_seq)
        1401

        In the feature's local coordinate system the first letter is 0, and the
        last letter is 1400 (the feature is 1401 letters long). The first exon
        runs from letter 0 to letter 153 (and has length 154):

        >>> exon0_seq = f.sub_features[0].extract(record.seq)
        >>> len(exon0_seq)
        154
        >>> str(exon0_seq) == str(feature_seq[0:154])
        True
        
        Now 154 is the first letter of the second exon, and 434 is its last
        letter. If you wanted to extract this exon from the feature's sequence
        you would use [154:435] like so:

        >>> exon1_seq = f.sub_features[1].extract(record.seq)
        >>> len(exon1_seq)
        281
        >>> str(exon1_seq) == str(feature_seq[154:435])
        True

        First let's try some position within the first exon (0 to 153),
        
        >>> f.get_parent_coord(0)
        3461
        >>> f.get_parent_coord(153)
        3614

        Now what happens with a local coordindate within the second exon? i.e.
        Letters 154 to 434 of the feature's sequence:

        >>> f.get_parent_coord(154)
        3697
        >>> f.get_parent_coord(434)
        3977

        To demonstrate this is working just as expected for single letters:
        
        >>> for x in range(len(f)):
        ...     assert feature_seq[x] == record.seq[f.get_parent_coord(x)]

        Note that the exon end points make slicing the parent directly hard!
        In the feature coordindates 154 is used as the end slice for the first
        exon, and the start slice point of the following exon. However, in
        the parent sequence you need 3615 and 3697 in these cases!

        >>> str(exon0_seq) == str(feature_seq[0:154]) == str(record.seq[3461:3615])
        True
        >>> str(exon1_seq) == str(feature_seq[154:435]) == str(record.seq[3697:3978])
        True

        To try to explain, notice:

        >>> f.get_parent_coord(154-1)+1
        3615
        >>> f.get_parent_coord(154)
        3697
        
        Sadley given the way that Python slicing works, this annoyance seems
        unavoidable. So, if you do want to get the sequence of just the exon
        the first two methods are preferable.

        Now a similar example of a join on the reverse strand,        

        >>> f = record.features[2]
        >>> print f.type, f.strand, len(f.sub_features)
        CDS -1 6
        >>> for sf in f.sub_features:
        ...    print sf.location, len(sf)
        [6616:6953] 337
        [7265:7351] 86
        [7463:7603] 140
        [7915:7998] 83
        [8086:8166] 80
        [8272:8368] 96
        >>> print f.location, len(f)
        [6616:8368] 822
        >>> feature_seq = f.extract(record.seq)
        >>> len(feature_seq)
        822
        >>> f.get_parent_coord(0)
        8367
        >>> f.get_parent_coord(95)
        8272
        >>> f.get_parent_coord(821)
        6616
        
        And a more indepth check - notice that we have to (reverse) complement
        each letter because the feature and the parent are on opposite strands:
        
        >>> from Bio.Seq import reverse_complement as rc
        >>> for x in range(len(f)):
        ...     assert feature_seq[x] == rc(record.seq[f.get_parent_coord(x)])
        """
        if local_coordinate < 0 or local_coordinate > len(self):
            #TODO - Can we handle this case nicely to mean upstream/downstream
            #of the feature? What if this would go outside the parent sequence
            #(or should be wrapped for a circular parent)?
            raise ValueError("Local co-ordinate should be within the feature")
        if self.sub_features:
            if self.strand == -1:
                #Reverse strand
                x = local_coordinate
                for f in self.sub_features[::-1]:
                    assert f.strand == -1
                    l = len(f)
                    if x < l:
                        return f.get_parent_coord(x)
                    #Not in this subfeature, try the next
                    x -= l
                assert False, "Problem getting co-ordinate via sub-features"
            else:
                #Forward strand (or mixed strand)
                x = local_coordinate
                for f in self.sub_features:
                    l = len(f)
                    if x < l:
                        return f.get_parent_coord(x)
                    #Not in this subfeature, try the next
                    x -= l
                assert False, "Problem getting co-ordinate via sub-features"
        elif self.strand == -1:
            return self.location.end.position + self.location.end.extension \
                   - local_coordinate - 1
        else:
            return self.location.start.position + local_coordinate
        
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
        gene -1 [1716:4347]
        tRNA -1 [1716:4347]
        
        Note that for a feature defined as a join of several subfeatures (e.g.
        the union of several exons) the gaps are not checked (e.g. introns).
        In this example, the tRNA location is defined in the GenBank file as
        complement(join(1717..1751,4311..4347)), so that position 1760 falls
        in the gap:

        >>> for f in record.features:
        ...     if 1760 in f:
        ...         print f.type, f.strand, f.location
        source 1 [0:154478]
        gene -1 [1716:4347]

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
        if self.sub_features:
            for f in self.sub_features:
                if value in f:
                    return True
            return False
        else:
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
    def __init__(self, start, end):
        """Specify the start and end of a sequence feature.

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

        """
        if isinstance(start, AbstractPosition):
            self._start = start
        else:
            self._start = ExactPosition(start)

        if isinstance(end, AbstractPosition):
            self._end = end
        else:
            self._end = ExactPosition(end)

    def __str__(self):
        """Returns a representation of the location (with python counting).

        For the simple case this uses the python splicing syntax, [122:150]
        (zero based counting) which GenBank would call 123..150 (one based
        counting).
        """
        return "[%s:%s]" % (self._start, self._end)

    def __repr__(self):
        """A string representation of the location for debugging."""
        return "%s(%s,%s)" \
               % (self.__class__.__name__, repr(self.start), repr(self.end))

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
        if value < self._start.position \
        or value >= self._end.position + self._end.extension:
            return False
        else:
            return True

    def _shift(self, offset):
        """Returns a copy of the location shifted by the offset (PRIVATE)."""
        return FeatureLocation(start = self._start._shift(offset),
                               end = self._end._shift(offset))

    start = property(fget= lambda self : self._start,
                 doc="Start location (possibly a fuzzy position, read only).")

    end = property(fget= lambda self : self._end,
                   doc="End location (possibly a fuzzy position, read only).")

    def _get_nofuzzy_start(self):
        #TODO - Do we still use the BetweenPosition class?
        if ((self._start == self._end) and isinstance(self._start,
             BetweenPosition)):
            return self._start.position
        else:
            return min(self._start.position,
                       self._start.position + self._start.extension)
    nofuzzy_start = property(fget=_get_nofuzzy_start,
        doc="""Start position (integer, approximated if fuzzy, read only).

        To get non-fuzzy attributes (ie. the position only) ask for
        'location.nofuzzy_start', 'location.nofuzzy_end'. These should return
        the largest range of the fuzzy position. So something like:
        (10.20)..(30.40) should return 10 for start, and 40 for end.
        """)

    def _get_nofuzzy_end(self):
        #TODO - Do we still use the BetweenPosition class?
        if ((self._start == self._end) and isinstance(self._start,
             BetweenPosition)):
            return self._end.position
        else:
            return max(self._end.position,
                       self._end.position + self._end.extension)
    nofuzzy_end = property(fget=_get_nofuzzy_end,
        doc="""End position (integer, approximated if fuzzy, read only).

        To get non-fuzzy attributes (ie. the position only) ask for
        'location.nofuzzy_start', 'location.nofuzzy_end'. These should return
        the largest range of the fuzzy position. So something like:
        (10.20)..(30.40) should return 10 for start, and 40 for end.
        """)

class AbstractPosition(object):
    """Abstract base class representing a position.
    """
    def __init__(self, position, extension):
        self.position = position
        self.extension = extension

    def __repr__(self):
        """String representation of the location for debugging."""
        return "%s(%s,%s)" % (self.__class__.__name__, \
                              repr(self.position), repr(self.extension))

    def __cmp__(self, other):
        """A simple comparison function for positions.

        This is very simple-minded and just compares the position attribute
        of the features; extensions are not considered at all. This could
        potentially be expanded to try to take advantage of extensions.
        """
        assert isinstance(other, AbstractPosition), \
          "We can only do comparisons between Biopython Position objects."

        return cmp(self.position, other.position)

    def _shift(self, offset):
        #We want this to maintain the subclass when called from a subclass
        return self.__class__(self.position + offset, self.extension)
            
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
