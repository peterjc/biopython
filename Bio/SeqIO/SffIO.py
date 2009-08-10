# Copyright 2009 by Jose Blanca (COMAV-UPV) and Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.SeqIO support for the Roche binary SFF file format.

You are expected to use this module via the Bio.SeqIO functions."""

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import struct

class _BinaryFile(object):
    '''Class for reading data from a binary file (PRIVATE).'''
    def __init__(self, fhandler, endianess='native'):
        '''It initializes the BinaryFile.
        
        keyword arguments:
        fhandler -- a file handler instance to the binary file.
        endianess -- native, big or little
        '''
        self._fh = fhandler
        self.bytes_read = None  #number of bytes read in the last time
        self._endianess = None
        self.endianess = endianess
    
    def _set_endianess(self, endianess):
        '''It sets the endianess that should be used to read.'''
        if endianess == 'native':
            self._endianess = '='
        elif endianess == 'big':
            self._endianess = '>'
        elif endianess == 'little':
            self._endianess = '<'
        else:
            raise ValueError('This endianess type is not supported')
    def _get_endianess(self):
        '''It returns the endianess of the file'''
        return self._endianess
    endianess = property(_get_endianess, _set_endianess)

    def seek(self, offset, whence=None):
        '''It moves to new file position.
        
        Argument offset is a byte count.  Optional argument whence defaults to
        0 (offset from start of file, offset should be >= 0); other values are
        1 (move relative to current position, positive or negative), and 2
        (move relative to end of file, usually negative, although many platforms
        allow seeking beyond the end of a file).
        keyword arguments:
        offset -- byte count
        whence -- (default 0)
        '''
        if whence is None:
            whence = 0
        self._fh.seek(offset, whence)
    def tell(self):
        '''It returns the current position in the file'''
        return self._fh.tell()

    def read(self, struct_def, byte_padding=None):
        '''It reads a chunk of a binary file.
    
        You have to provide the struct_definition (like '4c') and the byte_padding.
        If a byte_padding is given the number of bytes read will be a multiple of
        that number, adding the required pad at the end.
        It returns the values read as a tuple. It returns a tuple because the
        struct definition can have several things in it (like a char and an it
        ('ci'))
        It will alter the position of the file using the seek method.
        '''
        #Types to be used in the struct
        #Format    C Type                Python
        #  x      pad byte               no value   
        #  c      char                   string of length 1   
        #  b      signed char            integer   
        #  B      unsigned char          integer   
        #  ?      _Bool                  bool
        #  h      short                  integer   
        #  H      unsigned short         integer   
        #  i      int                    integer   
        #  I      unsigned int           integer or long   
        #  l      long                   integer   
        #  L      unsigned long          long   
        #  q      long long              long
        #  Q      unsigned long long     long
        #  f      float                  float   
        #  d      double                 float   
        #  s      char[]                 string   
        #  p      char[]                 string   
        #  P      void *                 long   
    
        offset = self._fh.tell()
        #if there is no explicit endianess given maybe there's a default one
        if self._endianess is not None and \
            struct_def[0] not in ('=', '@', '<', '>', '!'):
            struct_def = self._endianess + struct_def
        #we go to the place and read
        bytes_read = struct.calcsize(struct_def)
        chunk = self._fh.read(bytes_read)
        read = struct.unpack(struct_def, chunk)
    
        #if there is byte_padding the bytes_to_read should be a multiple of the
        #byte_padding
        if byte_padding is not None:
            pad = byte_padding
            bytes_read = ((bytes_read + pad - 1) // pad) * pad
            self._fh.seek(offset + bytes_read)
        self.bytes_read = bytes_read
        return read

    def read_record(self, record_definition, byte_padding=None):
        '''It reads a complete record (formed by several fields) in one go.
        
        The record should be defined as a list of tuples. Each tuple should have
        the field_name and the data type (and optionally the byte padding),
        like:
        record = [('1st','c'),('2ndr','c', 8)]
        It returns a dict with the field names as keys and the data read as values
        '''
        data = {}
        offset = self._fh.tell()
        bytes_read = 0
        for field in record_definition:
            if len(field) == 2: #no padding
                value = self.read(field[1])
            else: #padding
                value = self.read(field[1], field[2])
            bytes_read += self.bytes_read
            data[field[0]] = value
        #we have to sum up all the bytes reads for all fields
        #is there any padding?
        if byte_padding is not None:
            pad = byte_padding
            bytes_read = ((bytes_read + pad - 1) // pad) * pad
            self._fh.seek(offset + bytes_read)
        self.bytes_read = bytes_read
        return data

    def read_records(self, record_definition, num_records):
        '''It reads a list of records.
        
        All records should obey the given record_definition. num_records will be
        read.
        A generator with the read records will be returned.
        ''' 
        # pylint: disable-msg=W0612
        for i in range(num_records):    
            yield self.read_record(record_definition)

class _SffReader(object):
    'It reads sff files, returns data using dictionaries (PRIVATE).'
    #pylint: disable-msg=R0903
    def __init__(self, handle):
        'inits the SffReader class. It requires an sff file opened in binary mode.'
        if "b" not in handle.mode.lower() :
            raise ValueError("Requires the handle be opened in binary mode.")
        self._sff_bin = _BinaryFile(handle)
        #sff file are big endian
        self._sff_bin.endianess = 'big'
        self._number_of_reads = None
        self._header = {}   #The header data
        self._read_header()

    def _get_number_of_reads(self):
        'It returns the number of reads'
        return self._number_of_reads
    number_of_reads = property(_get_number_of_reads)

    @staticmethod
    def _check_magic_number(magic_number):
        'It checks that the magic number is the sff one'
        if magic_number != 779314790:
            raise RuntimeError('This file does not seems to be an sff file.')

    @staticmethod
    def _check_version(version):
        '''It checks that the version is supported, otherwise it raises an
        error.'''
        supported = ['\x00', '\x00', '\x00', '\x01']
        if version != supported:
            raise RuntimeError('SFF version not supported.')

    def _read_header(self):
        'It reads the info contained in the sff header'
        #the first part of the header is independent of the number of
        #the sequences
        #magic_number               I
        #version                    4c
        #index_offset               Q
        #index_length               I
        #number_of_reads            I
        #header_length              H
        #key_length                 H
        #number_of_flows_per_read   H
        #flowgram_format_code       B
        fmt = 'I4cQIIHHHB'
        self._sff_bin.seek(0)
        header = self._header
        (magic, ver0, ver1, ver2, ver3, header['index_offset'],
         header['index_length'], number_of_reads, header['header_length'],
         key_length, header['number_of_flows'], header['flowgram_code']) = \
                self._sff_bin.read(fmt)
        self._check_magic_number(magic)
        self._check_version([ver0, ver1, ver2, ver3])
        self._number_of_reads = number_of_reads

        #now we can read the rest of the header
        fmt = str(header['number_of_flows']) + 'c' + str(key_length) + 'c'
        data = self._sff_bin.read(fmt)
        header['flow_chars'] = data[:header['number_of_flows']]
        header['key_sequence'] = data[header['number_of_flows']:]

    def _read_sequence(self):
        'It reads one sequence from the sff file'
        bin = self._sff_bin
        read_offset = bin.tell()
        read = {} #here we'll store the info for this read
        #the sequence header format
        #read_header_length     H
        #name_length            H
        #number_of_bases        I
        #clip_qual_left         H
        #clip_qual_right        H
        #clip_adapter_left      H
        #clip_adapter_right     H
        header = self._header
        fmt = '>2HI4H'
        (head_length, name_length, num_bases, read['clip_qual_left'],
         read['clip_qual_right'], read['clip_adapter_left'],
         read['clip_adapter_right']) = bin.read(fmt)
        #now the name
        fmt = str(name_length) + 'c'
        read['name'] = ''.join(bin.read(fmt))
        #now the flowgram, bases and qualities
        if header['flowgram_code'] == 1:
            flow_type = 'H'
        else:
            raise RuntimeError('flowgram format code not supported')
        bin.seek(read_offset + head_length)
        num_bases = str(num_bases)
        flows_per_read = str(header['number_of_flows'])
        rec = (('flowgram_values', flows_per_read + flow_type),
               ('flow_index_per_base', num_bases + 'B'),
               ('bases', num_bases + 'c'),
               ('quality_scores', num_bases + 'B'))
        data = bin.read_record(rec, byte_padding=8)
        read['bases'] = ''.join(data['bases'])
        read['qualities'] = data['quality_scores']
        return read

    def sequences(self):
        'It yields the sequences from the sff file as SeqRecord'
        #where are the sequences?
        #after the header
        #in fact the documentation says that if there is an index (an
        #optional feature) the sequences should go after it, but that's not
        #the case in the test sequences I've got. In that case we should do
        #if not header['index_offset']:
        #    seqs_offset = header['index_length'] + header['index_offset']
        seqs_offset = self._header['header_length']
        self._sff_bin.seek(seqs_offset)
        n_reads = 0
        while n_reads < self.number_of_reads:
            read_dict = self._read_sequence()
            yield read_dict
            n_reads += 1

#This is a generator function!
def SffIterator(handle, alphabet = generic_dna, trim=False) :
    """Iterate over Roche 454 SFF reads (as SeqRecord objects).

    handle - input file, a Roche 454 file opened in binary mode.
    alphabet - optional alphabet, defaults to generic DNA.
    trim - should the sequences be trimmed?

    The resulting SeqRecord objects should match those from a paired
    FASTA and QUAL file converted from the SFF file using the Roche
    454 tool ssfinfo. i.e. The sequence will be mixed case, with the
    trim regions shown in lower case.
    """
    reader = _SffReader(handle) #Hack

    for read in reader.sequences():
        #Mimic the FASTA untrimmed output from the Roche tools
        seq = read["bases"]
        clip_qual_left = read["clip_qual_left"] - 1 #python slicing
        clip_qual_right = read["clip_qual_right"]
        assert 0 == read["clip_adapter_left"] and \
               0 == read["clip_adapter_right"], \
               "Cliping adators not supported yet" #Need test cases
        quals = list(read['qualities']) #Used lists in FASTQ and QUAL parsers
        if trim :
            seq = seq[clip_qual_left:clip_qual_right].upper()
            quals = quals[clip_qual_left:clip_qual_right]
        else :
            seq = seq[:clip_qual_left].lower() + \
                  seq[clip_qual_left:clip_qual_right].upper() + \
                  seq[clip_qual_right:].lower()
        record = SeqRecord(Seq(seq, alphabet),
                           id=read["name"],
                           name=read["name"])
        #for key in ["clip_qual_left", "clip_qual_right",
        #            "clip_adapter_left", "clip_adapter_right"] :
        #    record.annotations[key] = read[key]
        record.letter_annotations["phred_quality"] = quals
        #TODO - clipping
        #TODO - paired reads
        #Return the record and then continue...
        yield record

if __name__ == "__main__" :
    print "Running quick self test"
    filename = "../../Tests/Roche/E3MFGYR02_random_10_reads.sff"
    sff = list(SffIterator(open(filename,"rb")))
    sff_trim = list(SffIterator(open(filename,"rb"), trim=True))

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
