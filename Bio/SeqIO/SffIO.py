# Copyright 2009 by Peter Cock.  All rights reserved.
# Based on code contributed and copyright 2009 by Jose Blanca (COMAV-UPV).
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.SeqIO support for the binary Standard Flowgram Format (SFF) file format.

SFF was designed by 454 Life Sciences (Roche), the Whitehead Institute for
Biomedical Research and the Wellcome Trust Sanger Institute. You are expected
to use this module via the Bio.SeqIO functions under the format name "sff".

For a description of the file format, please see:
http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=formats

"""
#TODO - Can we parse the (optional) index?

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import struct

def _sff_file_header(handle) :
    """Read in an SFF file header (PRIVATE).

    Assumes the handle is at the start of the file, will read forwards
    though the header and leave the handle pointing at the first record.
    Returns a tuple of values from the header.
    """
    if hasattr(handle,"mode") and "U" in handle.mode.upper() :
        raise ValueError("SFF files must not be opened in universal new lines mode. "
                         "Default or binary mode is fine.")
    #file header (part one)
    #use big endiean encdoing   >
    #magic_number               I
    #version                    4B
    #index_offset               Q
    #index_length               I
    #number_of_reads            I
    #header_length              H
    #key_length                 H
    #number_of_flows_per_read   H
    #flowgram_format_code       B
    #[rest of file header depends on the number of flows and how many keys]
    fmt = '>I4BQIIHHHB'
    assert 31 == struct.calcsize(fmt)
    magic_number, ver0, ver1, ver2, ver3, index_offset, index_length, \
    number_of_reads, header_length, key_length, number_of_flows_per_read, \
    flowgram_format = struct.unpack(fmt, handle.read(31))
    if magic_number != 779314790 :
        raise ValueError("Wrong SFF magic number in header")
    if (ver0, ver1, ver2, ver3) != (0,0,0,1) :
        raise ValueError("Unsupported SFF version in header, %i.%i.%i.%i" \
                         % (ver0, ver1, ver2, ver3))
    if flowgram_format != 1 :
        raise ValueError("Flowgram format code %i not supported" \
                         % flowgram_format)
    if (index_offset!=0) ^ (index_length!=0) :
        raise ValueError("Index offset %i but index length %i" \
                         % (index_offset, index_length))
    flow_chars = handle.read(number_of_flows_per_read)
    key_sequence = handle.read(key_length)
    #According to the spec, the header_length field should be the total number
    #of bytes required by this set of header fields, and should be equal to
    #"31 + number_of_flows_per_read + key_length" rounded up to the next value
    #divisible by 8.
    assert header_length % 8 == 0
    padding = header_length - number_of_flows_per_read - key_length - 31
    assert 0 <= padding < 8, padding
    if chr(0)*padding != handle.read(padding) :
        raise ValueError("Post header %i byte padding region contained data" \
                         % padding)
    return header_length, index_offset, index_length, \
           number_of_reads, number_of_flows_per_read, \
           flow_chars, key_sequence

#This is a generator function!
def _sff_do_slow_index(handle) :
    """Generates an index by scanning though all the reads in an SFF file (PRIVATE).

    This is a slow but generic approach if we can't parse the provided index (if
    present).

    Will use the handle seek/tell functions.
    """
    handle.seek(0)
    header_length, index_offset, index_length, number_of_reads, \
    number_of_flows_per_read, flow_chars, key_sequence \
        = _sff_file_header(handle)
    #Now on to the reads...
    read_header_fmt = '>2HI4H'
    read_header_size = struct.calcsize(read_header_fmt)
    #NOTE - assuming flowgram_format==1, which means struct type H
    read_flow_fmt = ">%iH" % number_of_flows_per_read
    read_flow_size = struct.calcsize(read_flow_fmt)
    assert 1 == struct.calcsize(">B")
    assert 1 == struct.calcsize(">s")
    assert 1 == struct.calcsize(">c")
    assert read_header_size % 8 == 0 #Important for padding calc later!
    for read in range(number_of_reads) :
        record_offset = handle.tell()
        #assert record_offset%8 == 0 #Worth checking, but slow
        #First the fixed header
        read_header_length, name_length, seq_len, clip_qual_left, \
        clip_qual_right, clip_adapter_left, clip_adapter_right \
            = struct.unpack(read_header_fmt, handle.read(read_header_size))
        if read_header_length < 10 or read_header_length%8 != 0 :
            raise ValueError("Malformed read header, says length is %i" \
                             % read_header_length)
        #now the name and any padding (remainder of header)
        name = handle.read(name_length)
        padding = read_header_length - read_header_size - name_length
        if chr(0)*padding != handle.read(padding) :
            raise ValueError("Post name %i byte padding region contained data" \
                             % padding)
        assert record_offset + read_header_length == handle.tell()
        #now the flowgram values, flowgram index, bases and qualities
        size = read_flow_size + 3*seq_len
        handle.seek(size,1)
        #now any padding...
        padding = size%8
        if padding :
            padding = 8 - padding
            if chr(0)*padding != handle.read(padding) :
                raise ValueError("Post quality %i byte padding region contained data" \
                                 % padding)
        #print read, name, record_offset
        yield name, record_offset
    if handle.tell() % 8 != 0 :
        raise ValueError("After scanning reads, did not end on a multiple of 8")

def _sff_find_roche_index(handle) :
    """Locate any existing Roche style XML meta data and read index (PRIVATE).

    Makes a number of hard coded assumptions based on reverse engineered SFF
    files from Roche 454 machines.

    Returns a tuple of read count, SFF "index" offset and size, XML offset and
    size, and the actual read index offset and size."""    
    handle.seek(0)
    header_length, index_offset, index_length, number_of_reads, \
    number_of_flows_per_read, flow_chars, key_sequence \
        = _sff_file_header(handle)
    #print "Index offset %i, length %i, reads %i" \
    #      % (index_offset, index_length, number_of_reads)
    assert handle.tell() == header_length
    if not index_offset or not index_offset :
        raise ValueError("No index present in this SFF file")
    #Now jump to the header...
    handle.seek(index_offset)
    fmt = ">I4BLL"
    fmt_size = struct.calcsize(fmt)
    magic_number, ver0, ver1, ver2, ver3, xml_size, data_size \
                  = struct.unpack(fmt, handle.read(fmt_size))
    if magic_number != 778921588 :
        raise ValueError("Wrong magic number in SFF index header")
    if (ver0, ver1, ver2, ver3) != (49,46,48,48) :
        raise ValueError("Unsupported version in index header, %i.%i.%i.%i" \
                         % (ver0, ver1, ver2, ver3))
    if index_length != fmt_size + xml_size + data_size :
        raise ValueError("Problem understanding index header")
    if data_size != 20 * number_of_reads :
        raise ValueError("Expect index data block of %i bytes (20 bytes per read). "
                         "Got %i bytes" % (20 * number_of_reads, data_size))
    return number_of_reads, header_length, \
           index_offset, index_length, \
           index_offset + fmt_size, xml_size, \
           index_offset + fmt_size + xml_size, data_size

def _sff_read_roche_index_xml(handle) :
    """Reads any existing Roche style XML meta data in the SFF "index" (PRIVATE).

    Will use the handle seek/tell functions. Returns a string.
    """
    number_of_reads, header_length, index_offset, index_length, xml_offset, \
    xml_size, read_index_offset, read_index_size = _sff_find_roche_index(handle)
    handle.seek(xml_offset)
    return handle.read(xml_size)


#This is a generator function!
def _sff_read_roche_index(handle) :
    """Reads any existing Roche style read index provided in the SFF file (PRIVATE).

    Will use the handle seek/tell functions.

    Note: There are a number of hard coded assumptions here (e.g. the read names are
    14 characters), some of which could be relaxed given some suitable example files.
    """
    number_of_reads, header_length, index_offset, index_length, xml_offset, \
    xml_size, read_index_offset, read_index_size = _sff_find_roche_index(handle)
    #Now parse the read index...    
    handle.seek(read_index_offset)
    fmt = ">14s6B"
    assert 20 == struct.calcsize(fmt)
    for read in range(number_of_reads) :
        data = handle.read(20)
        name, x0, off3, off2, off1, off0, x255 = struct.unpack(fmt, data)
        if x0 != 0 :
            raise ValueError("Found %s instead of null at end of name for index entry %i,\n%s" \
                             % (repr(x0), read, repr(data)))
        if x255 != 255 :
            raise ValueError("Found %s instead of 0xff at end of index entry %i,\n%s" \
                             % (repr(x255), read, repr(data)))
        #TODO - Work out why struct doesn't do what I want with "L"
        offset = off0 + 255*off1 + 65025*off2 + 16581375*off3
        assert header_length <= offset <= index_offset, offset
        #print read, name, offset
        yield name, offset

def _sff_read_seq_record(handle, number_of_flows_per_read, alphabet,
                         trim=False) :
    """Parse the next read in the file, return data as a SeqRecord (PRIVATE)."""
    #Now on to the reads...
    #the read header format (fixed part):
    #read_header_length     H
    #name_length            H
    #seq_len                I
    #clip_qual_left         H
    #clip_qual_right        H
    #clip_adapter_left      H
    #clip_adapter_right     H
    #[rest of read header depends on the name length etc]
    read_header_fmt = '>2HI4H'
    read_header_size = struct.calcsize(read_header_fmt)
    read_flow_fmt = ">%iH" % number_of_flows_per_read
    read_flow_size = struct.calcsize(read_flow_fmt)

    read_header_length, name_length, seq_len, clip_qual_left, \
    clip_qual_right, clip_adapter_left, clip_adapter_right \
        = struct.unpack(read_header_fmt, handle.read(read_header_size))
    if clip_qual_left : clip_qual_left -= 1 #python counting
    if clip_adapter_left : clip_adapter_left -= 1 #python counting
    if read_header_length < 10 or read_header_length%8 != 0 :
        raise ValueError("Malformed read header, says length is %i" \
                         % read_header_length)
    #now the name and any padding (remainder of header)
    name = handle.read(name_length)
    padding = read_header_length - read_header_size - name_length
    if chr(0)*padding != handle.read(padding) :
        raise ValueError("Post name %i byte padding region contained data" \
                         % padding)
    #now the flowgram values, flowgram index, bases and qualities
    #NOTE - assuming flowgram_format==1, which means struct type H
    flow_values = struct.unpack(read_flow_fmt, handle.read(read_flow_size))
    temp_fmt = ">%iB" % seq_len # used for flow index and quals
    flow_index = struct.unpack(temp_fmt, handle.read(seq_len))
    seq = handle.read(seq_len)
    quals = list(struct.unpack(temp_fmt, handle.read(seq_len)))
    #now any padding...
    padding = (read_flow_size + seq_len*3)%8
    if padding :
        padding = 8 - padding
        if chr(0)*padding != handle.read(padding) :
            raise ValueError("Post quality %i byte padding region contained data" \
                             % padding)
    #Now build a SeqRecord
    if trim :
        seq = seq[clip_qual_left:clip_qual_right].upper()
        quals = quals[clip_qual_left:clip_qual_right]
    else :
        #This use of mixed case mimics the Roche SFF tool's FASTA output
        seq = seq[:clip_qual_left].lower() + \
              seq[clip_qual_left:clip_qual_right].upper() + \
              seq[clip_qual_right:].lower()
    record = SeqRecord(Seq(seq, alphabet),
                       id=name,
                       name=name,
                       description="")
    #Dirty trick to speed up this line:
    #record.letter_annotations["phred_quality"] = quals
    dict.__setitem__(record._per_letter_annotations,
                     "phred_quality", quals)
    #TODO - flow data
    #TODO - adaptor clipping
    #TODO - paired reads
    #Return the record and then continue...
    return record

#This is a generator function!
def SffIterator(handle, alphabet=generic_dna, trim=False) :
    """Iterate over Standard Flowgram Format (SFF) reads (as SeqRecord objects).

    handle - input file, an SFF file, e.g. from Roche 454 sequencing.
             This must NOT be opened in universal read lines mode!
    alphabet - optional alphabet, defaults to generic DNA.
    trim - should the sequences be trimmed?

    The resulting SeqRecord objects should match those from a paired
    FASTA and QUAL file converted from the SFF file using the Roche
    454 tool ssfinfo. i.e. The sequence will be mixed case, with the
    trim regions shown in lower case.
    """
    header_length, index_offset, index_length, number_of_reads, \
    number_of_flows_per_read, flow_chars, key_sequence \
        = _sff_file_header(handle)
    #Now on to the reads...
    #the read header format (fixed part):
    #read_header_length     H
    #name_length            H
    #seq_len                I
    #clip_qual_left         H
    #clip_qual_right        H
    #clip_adapter_left      H
    #clip_adapter_right     H
    #[rest of read header depends on the name length etc]
    read_header_fmt = '>2HI4H'
    read_header_size = struct.calcsize(read_header_fmt)
    read_flow_fmt = ">%iH" % number_of_flows_per_read
    read_flow_size = struct.calcsize(read_flow_fmt)
    assert 1 == struct.calcsize(">B")
    assert 1 == struct.calcsize(">s")
    assert 1 == struct.calcsize(">c")
    assert read_header_size % 8 == 0 #Important for padding calc later!
    #TODO - The spec allows for the index block to be in the middle
    #of the reads. We can check that if we keep track of our position
    #in the file...
    for read in range(number_of_reads) :
        yield _sff_read_seq_record(handle,
                                   number_of_flows_per_read,
                                   alphabet, trim)

#This is a generator function!
def _SffTrimIterator(handle, alphabet = generic_dna) :
    """Iterate over SFF reads (as SeqRecord objects) with trimming (PRIVATE)."""
    return SffIterator(handle, alphabet, trim=True)

if __name__ == "__main__" :
    print "Running quick self test"
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.sff"
    index1 = sorted(_sff_read_roche_index(open(filename)))
    index2 = sorted(_sff_read_roche_index(open(filename, "rB")))
    assert index1 == index2
    index2 = sorted(_sff_do_slow_index(open(filename)))
    assert index1 == index2
    index2 = sorted(_sff_do_slow_index(open(filename, "rB")))
    assert index1 == index2
    assert len(index1) == len(list(SffIterator(open(filename))))
    assert len(index1) == len(list(SffIterator(open(filename, "r"))))
    assert len(index1) == len(list(SffIterator(open(filename, "rB"))))
    from StringIO import StringIO
    assert len(index1) == len(list(SffIterator(StringIO(open(filename).read()))))
    assert len(index1) == len(list(SffIterator(StringIO(open(filename,"r").read()))))
    assert len(index1) == len(list(SffIterator(StringIO(open(filename,"rB").read()))))
                    
    sff = list(SffIterator(open(filename)))
    sff_trim = list(SffIterator(open(filename), trim=True))

    print _sff_read_roche_index_xml(open(filename))

    from Bio import SeqIO
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads_no_trim.fasta"
    fasta_no_trim = list(SeqIO.parse(open(filename,"rU"), "fasta"))
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads_no_trim.qual"
    qual_no_trim = list(SeqIO.parse(open(filename,"rU"), "qual"))

    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.fasta"
    fasta_trim = list(SeqIO.parse(open(filename,"rU"), "fasta"))
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.qual"
    qual_trim = list(SeqIO.parse(open(filename,"rU"), "qual"))

    for s, sT, f, q, fT, qT in zip(sff, sff_trim, fasta_no_trim, qual_no_trim, fasta_trim, qual_trim) :
        #print
        print s.id
        #print s.seq
        #print s.letter_annotations["phred_quality"]
        
        assert s.id == f.id == q.id
        assert str(s.seq) == str(f.seq)
        assert s.letter_annotations["phred_quality"] == q.letter_annotations["phred_quality"]

        assert s.id == sT.id == fT.id == qT.id
        assert str(sT.seq) == str(fT.seq)
        assert sT.letter_annotations["phred_quality"] == qT.letter_annotations["phred_quality"]

    print "Done"
