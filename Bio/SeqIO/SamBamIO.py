# Copyright 2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This module is for reading SAM (Sequence Alignment/Map) format files, and
# the closely related binary variant BAM format files, as SeqRecord objects.
# This module is expected to be used via the Bio.SeqIO API.

"""Bio.SeqIO support for the SAM and BAM file formats.

Note that you are expected to use this code via the Bio.SeqIO interface, as
shown below.

The SAM (Sequence Alignment/Map) format and its binary variant, BAM, are a
generic format for storing large nucleotide sequence alignments as generated
in second generation sequencing. They can also hold unaligned reads. The
Bio.SeqIO module in Biopython allows memory efficient access to the reads in
a SAM or BAM file, allowing for example easy conversion to FASTQ:

    >>> from Bio import SeqIO
    >>> SeqIO.convert("SamBam/ex1.sam", "sam", "temp_ex1.fastq", "fastq")
    3270

Or, from the BAM file instead:

    >>> from Bio import SeqIO
    >>> SeqIO.convert("SamBam/ex1.bam", "bam", "temp_ex1.fastq", "fastq")
    3270
    >>> import os
    >>> os.remove("temp_ex1.fastq")

Looking at the first read directly:

    >>> print SeqIO.parse("SamBam/ex1.sam", "sam").next().format("fastq")
    @EAS56_57:6:190:289:82/1
    CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA
    +
    <<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;
    <BLANKLINE>

    >>> print SeqIO.parse("SamBam/ex1.bam", "bam").next().format("fastq")
    @EAS56_57:6:190:289:82/1
    CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA
    +
    <<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;
    <BLANKLINE>

You can also get random access to the reads in a SAM or BAM file by their ID
(but there is no support for indexing by the mapping position):

    >>> from Bio import SeqIO
    >>> sam_reads = SeqIO.index("SamBam/ex1.sam", "sam")
    >>> print sam_reads["EAS56_57:6:190:289:82/2"].format("fastq")
    @EAS56_57:6:190:289:82/2
    AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC
    +
    <<<<<<;<<<<<<<<<<;<<;<<<<;8<6;9;;2;
    <BLANKLINE>

    >>> from Bio import SeqIO
    >>> bam_reads = SeqIO.index("SamBam/ex1.bam", "bam")
    >>> print bam_reads["EAS56_57:6:190:289:82/2"].format("fastq")
    @EAS56_57:6:190:289:82/2
    AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC
    +
    <<<<<<;<<<<<<<<<<;<<;<<<<;8<6;9;;2;
    <BLANKLINE>

For more indepth handling of SAM/BAM files, see the SAM tools project. Their
C API has a Python wrapper called pysam which is very helpful for more
complex task such as accessing columns of the alignment etc:
http://samtools.sourceforge.net
"""

import gzip
import struct

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq, _dna_complement_table
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import SANGER_SCORE_OFFSET

def _decode_flag(flag):
    paired = flag & 0x1
    if paired:
        well_mapped_pair = bool(flag & 0x2)
        mate_unmapped = bool(flag & 0x8)
        mate_strand = bool(flag & 0x20)
        pair1 = bool(flag & 0x40)
        pair2 = bool(flag & 0x80)
    else:
        well_mapped_pair = False
        assert not (flag & 0x2)
        mate_unmapped = False
        assert not flag & 0x8
        mate_rev_strand = False
        assert not flag & 0x20
        pair1 = False
        assert not flag & 0x40
        pair2 = False
        assert not flag & 0x80
        assert not (pair1 and pair2)
    unmapped = bool(flag & 0x4)
    rev_strand = bool(flag & 0x10)
    not_primary = bool(flag & 0x100)
    fails_qc = bool(flag & 0x200)
    duplicate = bool(flag & 0x400)
    return pair1, pair2

def SamIterator(handle, alphabet=generic_dna):
    #Precomputing FASTQ style encoding. Can we share this code with QualityIO?
    q_mapping = dict()
    for letter in range(0, 255):
        q_mapping[chr(letter)] = letter-SANGER_SCORE_OFFSET
    for line in handle:
        if line[0] == "@":
            #Ignore any optional header
            continue
        #Note that the tags field (and its tab character) are optional
        name, flag, ref_name, ref_pos, map_qual, cigar, \
        mate_reference_seq_name, mate_ref_pos, inferred_insert_size, \
        seq_string, quality_string, tags = (line.rstrip("\r\n")+"\t").split("\t",11)
        tags = tags.split("\t")
        
        qualities = [q_mapping[letter] for letter in quality_string]

        #TODO - How to handle any "." in the sequece?
        
        #SAM spec allows upper/lower case seq, while BAM uses a binary encoding
        #which ignores the case. The convention is to treat everything as
        #upper case (tolerate lower case in SAM but encourage upper case).
        yield _make_seq_record(name, seq_string.upper(), alphabet, qualities,
                               int(flag), int(map_qual))

def _make_seq_record(name, sequence, alphabet, qualities, flag, map_qual):
    identifier = name
    #Right now we don't need to decode everything, just basic pair information
    #pair1, pair2 = _decode_flag(flag)
    if flag & 0x040:
        identifier += "/1"
    elif flag & 0x080:
        identifier += "/2"
    if flag & 0x010:
        #Read was mapped onto the reverse strand, so stored in the SAM/BAM
        #file reverse complemented. I want the original read orientation...
        #Note for speed I assume the sequence is DNA
        return SeqRecord(Seq(sequence.translate(_dna_complement_table)[::-1],
                             alphabet),
                         id=identifier, name=name, description="",
                         annotations={"mapping_quality":map_qual},
                         letter_annotations={"phred_quality":qualities[::-1]})
    else:
        return SeqRecord(Seq(sequence, alphabet),
                         id=identifier, name=name, description="",
                         annotations={"mapping_quality":map_qual},
                         letter_annotations={"phred_quality":qualities})

def _bgzf_blocks(handle):
    r"""Parse the BGZF blocks (special GZIP header used in BAM files).
    
    The BAM file format uses a GZIP variant called BGZF. This consists of many
    GZIP blocks one after the other, with an extra subfield in each GZIP block
    header recording the block size. The GZIP format allows for such extra
    subfields. The BGZF standard is a subtype of GZIP using a specific
    compression (while GZIP in general allows several kinds) and a specific
    gzip flag (which indicates the presense of extra fields) where the block
    size is stored in an extra field under subfield identifier BC.
    
    Assumes the handle is at the start of the file and importantly has NOT
    been decompressed (e.g. with the gzip module).
    
    Returns block offset and size, plus the compressed data offset and size.
    Each block should be a self contained gzip entry with its own gzip header.

    >>> import gzip
    >>> handle = open("SamBam/ex1.bam", "rb")
    >>> for offsets in _bgzf_blocks(handle): print offsets
    (0, 18239, 18, 18214)
    (18239, 18223, 18257, 18198)
    (36462, 18017, 36480, 17992)
    (54479, 17342, 54497, 17317)
    (71821, 17715, 71839, 17690)
    (89536, 17728, 89554, 17703)
    (107264, 17292, 107282, 17267)
    (124556, 28, 124574, 3)

    Let's try and read the penultimate block, first check it starts with the
    expected header:

    >>> handle.seek(107264)
    >>> handle.read(4)
    '\x1f\x8b\x08\x04'

    Note this starts part way into a read. Now here is the clever bit :
    
    >>> block_offset = 107264
    >>> relative_offset = 101
    >>> handle.seek(block_offset)
    >>> h = gzip.GzipFile(fileobj=handle)
    >>> h.seek(relative_offset)
    >>> parts = _bam_file_read_header(h)
    >>> parts[0]
    'EAS139_19:7:85:262:751'
    >>> h.close()
    >>> handle.close()
    
    This is how the BAI (BAM index files) work. They give the BGZF/GZIP block
    offset (raw bytes in the compressed data) and the local offset within that
    block (bytes within the decompressed data).
    """
    while True:
        start_offset = handle.tell()
        magic = handle.read(4)
        if not magic:
            #End of file
            break
        if magic != "\x1f\x8b\x08\x04":
            raise ValueError(r"A BGZF (e.g. a BAM file) block should start with "
                             r"'\x1f\x8b\x08\x04' (decimal 31 139 8 4), not %s"
                             % repr(magic))
        gzip_mod_time = handle.read(4) #uint32_t
        gzip_extra_flags = handle.read(1) #uint8_t
        gzip_os = handle.read(1) #uint8_t
        extra_len = struct.unpack("<H", handle.read(2))[0] #uint16_t
        
        block_size = None
        x_len = 0
        while x_len < extra_len:
            subfield_id = handle.read(2)
            subfield_len = struct.unpack("<H", handle.read(2))[0] #uint16_t
            subfield_data = handle.read(subfield_len)
            x_len += subfield_len + 4
            if subfield_id == "BC":
                assert subfield_len == 2, "Wrong BC payload length"
                assert block_size is None, "Two BC subfields?"
                block_size = struct.unpack("<H", subfield_data)[0]+1 #uint16_t
        assert x_len == extra_len, (x_len, extra_len)
        #Now comes the compressed data, CRC, and length of uncompressed data.
        deflate_offset = handle.tell()
        deflate_size = block_size - extra_len - 19
        yield start_offset, block_size, deflate_offset, deflate_size
        #Note this will cope if the calling code moved the handle ;)
        handle.seek(start_offset + block_size)

from StringIO import StringIO #TODO - Use zlib instead of gzip+StringIO?
class _BgzfHandle(object):
    r"""Handle wrapper to read BGZF chunks (PRIVATE).
    
    >>> import gzip
    >>> h2 = gzip.open("SamBam/ex1.bam")
    >>> h2.seek(100)
    >>> data = h2.read(75000)
    >>> h2.tell()
    75100

    >>> handle = open("SamBam/ex1.bam", "rb")
    >>> h = _BgzfHandle(handle)
    >>> h.read(4)
    'BAM\x01'
    >>> h.tell()
    4
    >>> h.seek(4)
    >>> h.tell()
    4
    >>> h.seek(100)
    >>> h.tell()
    100
    >>> assert data == h.read(75000) #Should cover three gzip blocks
    >>> sorted(h._chunk_decompressed_offsets.iteritems())
    [(0, 0), (18239, 65536), (36462, 131072)]
    >>> assert len(data) == 75000, len(data)
    >>> h.tell()
    75100
    >>> h.seek(0)
    >>> h.tell()
    0
    >>> h.read(4)
    'BAM\x01'
    
    The clever bit of this is it allows random access to a BGFZ style GZIP file
    (i.e. random access to a BAM file), at the cost of decompressing each chunk
    into memory (typically under 20kb compressed data at a time).
    """
    def __init__(self, handle):
        self._handle = handle
        assert handle.tell()==0
        self._chunk_offset = 0
        self._chunk_decompressed_offsets = {0:0}
        next_chunk_offset, data = self._get_chunk(0)
        self._buffer = StringIO(data)
        #print "Decompressed data starts %s" % repr(data[:4])
        #print "Next offset %i" % self._next_chunk_offset

    def __iter__(self):
        raise NotImplementedError("Don't be silly, its binary not text")
    
    def readline(self):
        raise NotImplementedError("Don't be silly, its binary not text")
        
    def readlines(self):
        raise NotImplementedError("Don't be silly, that would be massive")
    
    def read(self, length):
        data = self._buffer.read(length)
        while len(data) < length:
            #Either end of file, or end of chunk.
            chunk_offsets = sorted(self._chunk_decompressed_offsets.keys())
            try:
                next_chunk_offset = chunk_offsets[chunk_offsets.index(self._chunk_offset)+1]
            except IndexError:
                assert self._chunk_offset == max(chunk_offsets)
                #Looks like end of file (and last chunk)
                return data
            try:
                next_next_chunk_offset, chunk_data = self._get_chunk(next_chunk_offset)
            except StopIteration:
                #Looks like end of file
                return data
            missing = length-len(data)
            self._buffer = StringIO(chunk_data)
            data += self._buffer.read(missing)
            self._relative_offset = 0          
            self._chunk_offset = next_chunk_offset
        return data
    
    def _get_chunk(self, start_offset):
        handle = self._handle
        handle.seek(start_offset)
        magic = handle.read(4)
        if not magic:
            #End of file
            #raise ValueError("End of file")
            raise StopIteration
        if magic != "\x1f\x8b\x08\x04":
            raise ValueError(r"A BGZF (e.g. a BAM file) block should start with "
                             r"'\x1f\x8b\x08\x04' (decimal 31 139 8 4), not %s"
                             % repr(magic))
        gzip_mod_time = handle.read(4) #uint32_t
        gzip_extra_flags = handle.read(1) #uint8_t
        gzip_os = handle.read(1) #uint8_t
        extra_len = struct.unpack("<H", handle.read(2))[0] #uint16_t
        
        block_size = None
        x_len = 0
        while x_len < extra_len:
            subfield_id = handle.read(2)
            subfield_len = struct.unpack("<H", handle.read(2))[0] #uint16_t
            subfield_data = handle.read(subfield_len)
            x_len += subfield_len + 4
            if subfield_id == "BC":
                assert subfield_len == 2, "Wrong BC payload length"
                assert block_size is None, "Two BC subfields?"
                block_size = struct.unpack("<H", subfield_data)[0]+1 #uint16_t
        assert x_len == extra_len, (x_len, extra_len)
        #Now comes the compressed data, CRC, and length of uncompressed data.
        deflate_offset = handle.tell()
        deflate_size = block_size - extra_len - 19
        #TODO - Should be able to use zlib instead of gzip...
        handle.seek(start_offset)
        data = gzip.GzipFile(fileobj=StringIO(handle.read(block_size))).read()
        
        self._chunk_decompressed_offsets[start_offset + block_size] \
            = self._chunk_decompressed_offsets[start_offset] + len(data)
        return start_offset + block_size, data
    
    def seek(self, offset):
        if self._chunk_decompressed_offsets[self._chunk_offset] <= offset \
        and offset < self._chunk_decompressed_offsets[self._chunk_offset] \
        + len(self._buffer.getvalue()):
            #Within current chunk
            self._buffer.seek(offset - self._chunk_decompressed_offsets[self._chunk_offset])
            assert offset == self.tell()
            return
        #Assume we've loaded enough chunks
        chunks = sorted(self._chunk_decompressed_offsets.iteritems())
        while chunks:
            chunk_offset, chunk_decompressed_offset = chunks.pop()
            if chunk_decompressed_offset <= offset:
                #Found the right chunk
                self._seek_chunk(chunk_offset, offset-chunk_decompressed_offset)
                break
        assert offset == self.tell()
    
    def tell(self):
        return self._chunk_decompressed_offsets[self._chunk_offset] \
                   + self._buffer.tell()
        #raise NotImplementedError("Use the tell_chunk method")
    
    def _seek_chunk(self, chunk_offset, relative_offset):
        if chunk_offset < 0:
            raise ValueError("Chunk offsets must be non-negative")
        if chunk_offset not in self._chunk_decompressed_offsets:
            raise ValueError("Haven't seen this chunk yet %i, only %s" % \
                             (chunk_offset,
                              repr(sorted(self._chunk_decompressed_offsets.keys()))))
        if relative_offset < 0:
            raise ValueError("Relative offsets must be non-negative")
        if chunk_offset != self._chunk_offset:
            try:
                next_chunk_offset, data = self._get_chunk(chunk_offset)
            except StopIteration:
                #End of file
                data = ""
            if len(data) < relative_offset:
                raise ValueError("Relative offset bigger than chunk")
            self._buffer = StringIO(data)
            self._chunk_offset = chunk_offset
        self._buffer.seek(relative_offset)
    
    def _tell_chunk(self):
        return self._chunk_offset, self._buffer.tell()
        

def _bgzf_index(handle):
    """Index a BAM file using the BGZF block structure.

    >>> handle = open("SamBam/ex1.bam", "rb")
    >>> data = list(_bgzf_index(handle))
    >>> data[0]
    ('EAS56_57:6:190:289:82/1', 38)
    >>> data[-1]
    ('EAS114_26:7:37:79:581/1', 456475)
    
    
    >>> h = _BgzfHandle(open("SamBam/ex1.bam", "rb"))
    >>> data = h.read(100)
    >>> while data: data = h.read(100)
    >>> h.tell()
    456614

    >>> h.seek(38)
    >>> h.tell()
    38
    >>> h.seek(456475)
    >>> h.tell()
    456475
    >>> parts = _bam_file_read_header(h)
    >>> parts[0]
    'EAS114_26:7:37:79:581'

    >>> sorted(h._chunk_decompressed_offsets.iteritems())
    [(0, 0), (18239, 65536), (36462, 131072), (54479, 196608), (71821, 262144), (89536, 327680), (107264, 393216), (124556, 456614), (124584, 456614)]
    >>> h._seek_chunk(124584, 0)
    >>> h._buffer.tell()
    0
    >>> h.tell()
    456614
    >>> for chunk, offset in sorted(h._chunk_decompressed_offsets.iteritems()):
    ...    print chunk, offset
    ...    h._seek_chunk(chunk, 0)
    ...    assert offset == h.tell()
    0 0
    18239 65536
    36462 131072
    54479 196608
    71821 262144
    89536 327680
    107264 393216
    124556 456614
    124584 456614

    """
    h = _BgzfHandle(handle)
    header, ref_count = _bam_file_header(h)
    #Skip any reference information
    for i in range(ref_count):
        ref_name, ref_len = _bam_file_reference(h)
    
    while True:
        try:
            key, start_offset, end_offset, ref_id, ref_pos, bin, \
                map_qual, cigar_len, flag, read_len, mate_ref_id, \
                mate_ref_pos = _bam_file_read_header(h)
        except StopIteration:
            #End of the reads
            break

        h.seek(start_offset)
        assert h.tell() == start_offset
        parts = _bam_file_read_header(h)
        assert parts[0] == key

        if flag & 0x40:
            key += "/1"
        elif flag & 0x80:
            key += "/2"
        yield key, start_offset
        h.seek(end_offset)


def _bam_file_header(handle):
    """Read in a BAM file header (PRIVATE).

    Assumings the handle is at the start of the file and has already been
    decompressed (e.g. with the gzip module).

    >>> handle = gzip.open("SamBam/ex1.bam")
    >>> header, num_refs = _bam_file_header(handle)
    >>> header
    ''
    >>> num_refs
    2
    """
    magic = handle.read(4)
    if magic != "BAM\1":
        raise ValueError("After decompression BAM files should start "
                         "with 'BAM\1', not %s" % repr(magic))
    assert 4 == struct.calcsize("<i")
    data = handle.read(4)
    #raise ValueError("Got %s" % repr(data))
    header_length = struct.unpack("<i", data)[0]
    header = handle.read(header_length).rstrip("\0")
    data = handle.read(4)
    #raise ValueError("Got %s" % repr(data))
    num_refs = struct.unpack("<i", data)[0]
    return header, num_refs

def _bam_file_reference(handle):
    """Parse next reference in a BAM file (PRIVATE).

    Assumings the handle is just after the header and has already been                                             
    decompressed (e.g. with the gzip module)

    >>> handle = gzip.open("SamBam/ex1.bam")
    >>> header, num_refs = _bam_file_header(handle)
    >>> for i in range(num_refs):
    ...     print _bam_file_reference(handle)
    ('chr1', 1575)
    ('chr2', 1584)

    """
    ref_name_len = struct.unpack("<i", handle.read(4))[0]
    ref_name = handle.read(ref_name_len).rstrip("\0")
    ref_len = struct.unpack("<i", handle.read(4))[0]
    return ref_name, ref_len

def _build_decoder():
    #TODO - Extend this to cover the ambiguous bases
    answer = {}
    for i,first in enumerate("ACGTN"):
        if first=="N":
            #16+32+64+128 = 240
            value = 240
        else:
            value = 2**(i+4)
        answer[value] = first+"-"
        for j,second in enumerate("ACGTN"):
            if second=="N":
                #1+2+4+8 = 15
                answer[value + 15] = first+second
            else:
                answer[value + 2**j] = first+second
    assert answer[79] == "GN"
    assert answer[240] == "N-"
    assert answer[241] == "NA"
    return answer
_decode_dibase_byte = _build_decoder()

def _bam_file_read_header(handle):
    """Parse the header of the next read in a BAM file (PRIVATE).

    >>> handle = gzip.open("SamBam/ex1.bam")
    >>> header, num_refs = _bam_file_header(handle)
    >>> for i in range(num_refs):
    ...     print _bam_file_reference(handle)
    ('chr1', 1575)
    ('chr2', 1584)
    >>> print _bam_file_read_header(handle)[:3]
    ('EAS56_57:6:190:289:82', 38, 153)
    >>> handle.seek(153)
    >>> print _bam_file_read_header(handle)[:3]
    ('EAS56_57:6:190:289:82', 153, 292)
    >>> handle.seek(153)
    >>> print _bam_file_read_header(handle)[:3]
    ('EAS56_57:6:190:289:82', 153, 292)

    Returns a tuple of the read name, start offset, end offset, etc.

    The offset information is used for indexing - the start offset is to find
    the record on demand, the end offset is used when building the index to
    skip over the rest of the read.
    
    The end offset is also used to determine the end of the tags section.
    """
    #TODO - Check BBH really works for the bin_mq_ml field, defined by
    # bin_mq_nl = bin<<16|mapQual<<8|read_name_len (including NULL)
    start_offset = handle.tell()
    
    fmt = "<iiiBBHHHiiii"
    assert 36 == struct.calcsize(fmt)
    data = handle.read(36)
    if not data:
        raise StopIteration
    if len(data) < 26:
        raise ValueError("Premature end of file")
    #raise ValueError("Data %s = %s" % (repr(data), repr(struct.unpack(fmt, data))))

    block_size, ref_id, ref_pos, read_name_len, map_qual, bin, \
    cigar_len, flag, read_len, mate_ref_id, mate_ref_pos, \
    inferred_insert_size = struct.unpack(fmt, data)

    if read_name_len > 50 or read_name_len <= 0:
        raise ValueError("A read name length of %i probably means the "
                         "parser is out of sync somehow. Read starts:\n%s"
                         "\nand the read name would be %s etc." \
                         % (read_name_len, repr(data), repr(handle.read(25))))
    if read_len > 5000 or read_len <= 0:
        raise ValueError("A read length of %i probably means the parser is out "
                         "of sync somehow. Read starts:\n%s\nand the read name "
                         "would be %s." \
                         % (read_len, repr(data), repr(handle.read(read_name_len))))

    read_name = handle.read(read_name_len).rstrip("\0")
    end_offset = start_offset + block_size + 4
    return read_name, start_offset, end_offset, ref_id, ref_pos, bin, \
           map_qual, cigar_len, flag, read_len, mate_ref_id, mate_ref_pos


def _bam_file_read(handle):
    """Parse the next read in a BAM file (PRIVATE).

    >>> handle = gzip.open("SamBam/ex1.bam")
    >>> header, num_refs = _bam_file_header(handle)
    >>> for i in range(num_refs):
    ...     print _bam_file_reference(handle)
    ('chr1', 1575)
    ('chr2', 1584)
    >>> print _bam_file_read(handle)[:2]
    ('EAS56_57:6:190:289:82', 'CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA')
    >>> print _bam_file_read(handle)[:2]
    ('EAS56_57:6:190:289:82', 'AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC')
    >>> print _bam_file_read(handle)[:2]
    ('EAS51_64:3:190:727:308', 'GGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGG')

    """
    read_name, start_offset, end_offset, ref_id, ref_pos, \
        bin, map_qual, cigar_len, flag, read_len, mate_ref_id, \
        mate_ref_pos = _bam_file_read_header(handle)
    
    fmt = "<%iI" % cigar_len
    cigar = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))
    #raise ValueError("flag %i, cigar length %i, %s" % (flag, cigar_len, repr(cigar)))

    fmt = "<%iB" % ((read_len+1)/2) # round up to make it even
    try:
        binary_seq = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))
    except Exception, err:
        raise ValueError("Expected read name %s of length %i\n%s" \
                         % (repr(read_name), read_len, err))
    seq_string = "".join(_decode_dibase_byte[byte] for byte in binary_seq)
    seq_string = seq_string[:read_len] # remove extra value if odd

    fmt = "<%iB" % read_len
    quals = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))

    #TODO - Parse the tags
    #For now, just skip the tags
    handle.seek(end_offset)

    return read_name, seq_string, quals, flag, map_qual

def BamIterator(handle, alphabet=generic_dna):
    h = gzip.GzipFile(fileobj=handle)
    header, ref_count = _bam_file_header(h)
    #Skip any reference information
    for i in range(ref_count):
        ref_name, ref_len = _bam_file_reference(h)
    #Loop over the reads
    while True:
        name, seq_string, qualities, flag, map_qual = _bam_file_read(h)
        yield _make_seq_record(name, seq_string, alphabet, qualities, flag,
                               map_qual)
    raise StopIteration

def _test():
    """Run the Bio.SeqIO.BamSamIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..", "..", "Tests")):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..", "..", "Tests"))
        assert os.path.isfile("SamBam/ex1.sam")
        assert os.path.isfile("SamBam/ex1.bam")
        doctest.testmod(verbose=0)
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
    elif os.path.isdir("Tests"):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir("Tests")
        assert os.path.isfile("SamBam/ex1.sam")
        assert os.path.isfile("SamBam/ex1.bam")
        doctest.testmod(verbose=0)
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
        
if __name__ == "__main__":
    _test()
