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
    return number_of_reads, number_of_flows_per_read

#This is a generator function!
def SffIterator(handle, alphabet = generic_dna, trim=False) :
    """Iterate over Standard Flowgram Format (SFF) reads (as SeqRecord objects).

    handle - input file, an SFF file, e.g. from Roche 454 sequencing.
    alphabet - optional alphabet, defaults to generic DNA.
    trim - should the sequences be trimmed?

    The resulting SeqRecord objects should match those from a paired
    FASTA and QUAL file converted from the SFF file using the Roche
    454 tool ssfinfo. i.e. The sequence will be mixed case, with the
    trim regions shown in lower case.
    """
    number_of_reads, number_of_flows_per_read = _sff_file_header(handle)
    #Now on to the reads...
    #the read header format (fixed part):
    #read_header_length     H
    #name_length            H
    #number_of_bases        I
    #clip_qual_left         H
    #clip_qual_right        H
    #clip_adapter_left      H
    #clip_adapter_right     H
    #[rest of read header depends on the name length etc]
    read_header_fmt = '>2HI4H'
    read_header_size = struct.calcsize(read_header_fmt)
    H_size = struct.calcsize(">H")
    assert 1 == struct.calcsize(">B")
    assert 1 == struct.calcsize(">s")
    assert 1 == struct.calcsize(">c")
    assert read_header_size % 8 == 0
    for read in range(number_of_reads) :
        #First the fixed header
        read_header_length, name_length, number_of_bases, clip_qual_left, \
        clip_qual_right, clip_adapter_left, clip_adapter_right \
            = struct.unpack(read_header_fmt, handle.read(read_header_size))
        if clip_qual_left : clip_qual_left -= 1 #python counting
        if clip_adapter_left : clip_adapter_left -= 1 #python counting
        if read_header_length < 10 :
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
        flow_values = struct.unpack(">%iH" % number_of_flows_per_read,
                                    handle.read(number_of_flows_per_read*H_size))
        flow_index = struct.unpack(">%iB" % number_of_bases,
                                   handle.read(number_of_bases))
        seq = handle.read(number_of_bases)
        quals = list(struct.unpack(">%iB" % number_of_bases,
                                   handle.read(number_of_bases)))
        #now any padding...
        padding = (number_of_flows_per_read*H_size + number_of_bases*3)%8
        if padding :
            padding = 8 - padding
            if chr(0)*padding != handle.read(padding) :
                raise ValueError("Post quality %i byte padding region contained data" \
                                 % padding)
        #Yield this read as a SeqRecord :)
        if trim :
            seq = seq[clip_qual_left:clip_qual_right].upper()
            quals = quals[clip_qual_left:clip_qual_right]
        else :
            seq = seq[:clip_qual_left].lower() + \
                  seq[clip_qual_left:clip_qual_right].upper() + \
                  seq[clip_qual_right:].lower()
        record = SeqRecord(Seq(seq, alphabet),
                           id=name,
                           name=name,
                           description="")
        record.letter_annotations["phred_quality"] = quals
        #TODO - flow data
        #TODO - adaptor clipping
        #TODO - paired reads
        #Return the record and then continue...
        yield record

#This is a generator function!
def _SffTrimIterator(handle, alphabet = generic_dna) :
    """Iterate over SFF reads (as SeqRecord objects) with trimming (PRIVATE)."""
    return SffIterator(handle, alphabet, trim=True)

if __name__ == "__main__" :
    print "Running quick self test"
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.sff"
    sff = list(SffIterator(open(filename)))
    sff_trim = list(SffIterator(open(filename), trim=True))

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
        print
        print s.id
        print s.seq
        print s.letter_annotations["phred_quality"]
        
        assert s.id == f.id == q.id
        assert str(s.seq) == str(f.seq)
        assert s.letter_annotations["phred_quality"] == q.letter_annotations["phred_quality"]

        assert s.id == sT.id == fT.id == qT.id
        assert str(sT.seq) == str(fT.seq)
        assert sT.letter_annotations["phred_quality"] == qT.letter_annotations["phred_quality"]
    print "Done"

