# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Dictionary like indexing of sequence files (PRIVATE).

You are not expected to access this module, or any of its code, directly. This
is all handled internally by the Bio.SeqIO.indexed_dict(...) function which is
the public interface for this.

The basic idea is that we scan over a sequence file, looking for new record
markers. We then try and extract the string that Bio.SeqIO.parse/read would
use as the record id, ideally without actually parsing the full record. We
then use a normal Python dictionary to record the file offset for the record
start against the record id.

Note that this means full parsing is on demand, so any invalid or problem record
may not trigger an exception until it is accessed. This is by design.

This means our dictionary like objects have in memory ALL the keys (all the
record identifiers), which shouldn't be a problem even with second generation
sequencing. If this is an issue later on, storing the keys and offsets in a
temp lookup file might be one idea.
"""

import re
from Bio import SeqIO

class _IndexedSeqFileDict(object) :
    """Read only dictionary interface to a sequential sequence file.

    Keeps the keys in memory, reads the file to access entries as
    SeqRecord objects using Bio.SeqIO for parsing them. Note - as
    with the Bio.SeqIO.to_dict() function, duplicate keys (record
    identifiers) are not allowed. If this happens, a ValueError
    exception is raised.
    """
    def __init__(self, filename, alphabet, mode="rU") :
        self._handle = open(filename, mode)
        self._index = dict()
        self._alphabet = alphabet
        self._format = ""
        #Now scan it in a subclassed method, and set the format!

    def _record_key(self, key, seek_position) :
        if key in self._index :
            raise ValueError("Duplicate key '%s'" % key)
        else :
            self._index[key] = seek_position

    def keys(self) :
        return self._index.keys()

    def __len__(self) :
        return len(self._index)

    def __contains__(self, value) :
        return value in self._index

    def __getitem__(self, key) :
        #For non-trivial file formats this must be over-ridden in the subclass
        handle = self._handle
        handle.seek(self._index[key])
        record = SeqIO.parse(handle, self._format, self._alphabet).next()
        assert record.id == key, "Requested key %s, found record.id %s" % (key, record.id)
        return record

    def get(self, k, d=None) :
        try :
            return self[k]
        except KeyError :
            return d

####################
# Special indexers #
####################

class SffDict(_IndexedSeqFileDict) :
    """Indexed dictionary like access to a Standard Flowgram Format (SFF) file."""
    def __init__(self, filename, alphabet) :
        _IndexedSeqFileDict.__init__(self, filename, alphabet, "r")
        handle = self._handle
        header_length, index_offset, index_length, number_of_reads, \
        self._flows_per_read = SeqIO.SffIO._sff_file_header(handle)
        if index_offset and index_length:
            #There is an index provided, try this the fast way:
            try :
                for name, offset in SeqIO.SffIO._sff_read_roche_index(handle) :
                    self._record_key(name, offset)
                assert len(self) == number_of_reads, \
                       "Indexed %i records, expected %i" \
                       % (len(self), number_of_reads)
                return
            except ValueError, err :
                import warnings
                warnings.warn("Could not parse the SFF index: %s" % err)
                self._index = {} #reset in case partially populated
                handle.seek(0)
        else :
            #TODO - Remove this debug warning?
            import warnings
            warnings.warn("No SFF index, doing it the slow way")
        #Fall back on the slow way!
        for name, offset in SeqIO.SffIO._sff_do_slow_index(handle) :
            #print "%s -> %i" % (name, offset)
            self._record_key(name, offset)
        assert len(self) == number_of_reads, \
               "Indexed %i records, expected %i" % (len(self), number_of_reads)

    def __getitem__(self, key) :
        handle = self._handle
        handle.seek(self._index[key])
        record = SeqIO.SffIO._sff_read_seq_record(handle,
                                                  self._flows_per_read,
                                                  self._alphabet)
        assert record.id == key
        return record

###################
# Simple indexers #
###################

class _SequentialSeqFileDict(_IndexedSeqFileDict) :
    """Subclass for easy cases (PRIVATE)."""
    def __init__(self, filename, alphabet, format, marker) :
        _IndexedSeqFileDict.__init__(self, filename, alphabet)
        self._format = format
        handle = self._handle
        marker_re = re.compile("^%s" % marker)
        marker_offset = len(marker)
        while True :
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            if marker_re.match(line) :
                #Here we can assume the record.id is the first word after the
                #marker. This is generally fine... but not for GenBank, EMBL, Swiss
                self._record_key(line[marker_offset:].strip().split(None,1)[0], offset)

class FastaDict(_SequentialSeqFileDict) :
    """Indexed dictionary like access to a FASTA file."""
    def __init__(self, filename, alphabet) :
        _SequentialSeqFileDict.__init__(self, filename, alphabet, "fasta", ">")

class QualDict(_SequentialSeqFileDict) :
    """Indexed dictionary like access to a QUAL file."""
    def __init__(self, filename, alphabet) :
        _SequentialSeqFileDict.__init__(self, filename, alphabet, "qual", ">")

class PirDict(_SequentialSeqFileDict) :
    """Indexed dictionary like access to a PIR/NBRF file."""
    def __init__(self, filename, alphabet) :
        _SequentialSeqFileDict.__init__(self, filename, alphabet, "pir", ">..;")

class PhdDict(_SequentialSeqFileDict) :
    """Indexed dictionary like access to a PHD (PHRED) file."""
    def __init__(self, filename, alphabet) :
        _SequentialSeqFileDict.__init__(self, filename, alphabet, "phd",
                                        "BEGIN_SEQUENCE")

class AceDict(_SequentialSeqFileDict) :
    """Indexed dictionary like access to an ACE file."""
    def __init__(self, filename, alphabet) :
        _SequentialSeqFileDict.__init__(self, filename, alphabet, "ace", "CO ")


#######################################
# Fiddly indexers: GenBank, EMBL, ... #
#######################################

class GenBankDict(_IndexedSeqFileDict) :
    """Indexed dictionary like access to a GenBank file."""
    def __init__(self, filename, alphabet) :
        _IndexedSeqFileDict.__init__(self, filename, alphabet)
        self._format = "genbank"
        handle = self._handle
        marker_re = re.compile("^LOCUS ")
        while True :
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            if marker_re.match(line) :
                #We cannot assume the record.id is the first word after LOCUS,
                #normally the first entry on the VERSION or ACCESSION line is used.
                key = None
                done = False
                while not done :
                    line = handle.readline()
                    if line.startswith("ACCESSION ") :
                        key = line.rstrip().split()[1]
                    elif line.startswith("VERSION ") :
                        version_id = line.rstrip().split()[1]
                        if version_id.count(".")==1 and version_id.split(".")[1].isdigit() :
                            #This should mimics the GenBank parser...
                            key = version_id
                            done = True
                            break
                    elif line.startswith("FEATURES ") \
                    or line.startswith("ORIGIN ") \
                    or line.startswith("//") \
                    or marker_re.match(line) \
                    or not line:
                        done = True
                        break
                if not key :
                    raise ValueError("Did not find ACCESSION/VERSION lines")
                self._record_key(key, offset)

class EmblDict(_IndexedSeqFileDict) :
    """Indexed dictionary like access to an EMBL file."""
    def __init__(self, filename, alphabet) :
        _IndexedSeqFileDict.__init__(self, filename, alphabet)
        self._format = "embl"
        handle = self._handle
        marker_re = re.compile("^ID ")
        while True :
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            if marker_re.match(line) :
                #We cannot assume the record.id is the first word after ID,
                #normally the SV line is used.
                parts = line[3:].rstrip().split(";")
                if parts[1].strip().startswith("SV ") :
                    #The SV bit gives the version
                    key = "%s.%s" % (parts[0].strip(),parts[1].strip().split()[1])
                else :
                    key = parts[0].strip()
                while True :
                    line = handle.readline()
                    if line.startswith("SV ") :
                        key = line.rstrip().split()[1]
                        break
                    elif line.startswith("FH ") \
                    or line.startswith("FT ") \
                    or line.startswith("SQ ") \
                    or line.startswith("//") \
                    or marker_re.match(line) \
                    or not line :
                        break
                self._record_key(key, offset)

class SwissDict(_IndexedSeqFileDict) :
    """Indexed dictionary like access to a SwissProt file."""
    def __init__(self, filename, alphabet) :
        _IndexedSeqFileDict.__init__(self, filename, alphabet)
        self._format = "swiss"
        handle = self._handle
        marker_re = re.compile("^ID ")
        while True :
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            if marker_re.match(line) :
                #We cannot assume the record.id is the first word after ID,
                #normally the following AC line is used.
                line = handle.readline()
                assert line.startswith("AC ")
                key = line[3:].strip().split(";")[0].strip()
                self._record_key(key, offset)

class IntelliGeneticsDict(_IndexedSeqFileDict) :
    """Indexed dictionary like access to a IntelliGenetics file."""
    def __init__(self, filename, alphabet) :
        _IndexedSeqFileDict.__init__(self, filename, alphabet)
        self._format = "ig"
        handle = self._handle
        marker_re = re.compile("^;")
        while True :
            offset = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            if marker_re.match(line) :
                #Now look for the first line which doesn't start ";"
                while True :
                    line = handle.readline()
                    if not line :
                        raise ValueError("Premature end of file?")
                    if line[0] != ";" and line.strip() :
                        key = line.split()[0]
                        self._record_key(key, offset)
                        break

##########################
# Now the FASTQ indexers #
##########################
         
class _FastqSeqFileDict(_IndexedSeqFileDict) :
    """Subclass for easy cases (PRIVATE).

    With FASTQ the records all start with a "@" line, but so too can some
    quality lines. For an initial implementation, this will only deal with
    non-line-wrapped FASTQ files (four lines per record).
    """
    def __init__(self, filename, alphabet, fastq_format) :
        _IndexedSeqFileDict.__init__(self, filename, alphabet)
        self._format = fastq_format
        handle = self._handle
        while True :
            pos = handle.tell()
            line = handle.readline()
            if not line : break #End of file
            if line[0] != "@" :
                raise ValueError("Problem with FASTQ @ line:\n%s" % repr(line))
            #This record seems OK (so far)
            self._record_key(line[1:].rstrip().split(None,1)[0],pos)
            handle.readline() #seq
            line = handle.readline()
            if line[0] != "+" :
                raise ValueError("Problem with FASTQ + line:\n%s" % repr(line))
            handle.readline() #qual

class FastqSangerDict(_FastqSeqFileDict) :
    """Indexed dictionary like access to a standard Sanger FASTQ file."""
    def __init__(self, filename, alphabet) :
        _FastqSeqFileDict.__init__(self, filename, alphabet, "fastq-sanger")

class FastqSolexaDict(_FastqSeqFileDict) :
    """Indexed dictionary like access to a Solexa (or early Illumina) FASTQ file."""
    def __init__(self, filename, alphabet) :
        _FastqSeqFileDict.__init__(self, filename, alphabet, "fastq-solexa")

class FastqIlluminaDict(_FastqSeqFileDict) :
    """Indexed dictionary like access to a Illumina 1.3+ FASTQ file."""
    def __init__(self, filename, alphabet) :
        _FastqSeqFileDict.__init__(self, filename, alphabet, "fastq-illumina")

###############################################################################

_FormatToIndexedDict = {"ace" : AceDict,
                        "embl" : EmblDict,
                        "fasta" : FastaDict,
                        "fastq" : FastqSangerDict,
                        "fastq-sanger" : FastqSangerDict, #alias of the above
                        "fastq-solexa" : FastqSolexaDict,
                        "fastq-illumina" : FastqIlluminaDict,
                        "genbank" : GenBankDict,
                        "gb" : GenBankDict, #alias of the above
                        "ig" : IntelliGeneticsDict,
                        "phd" : PhdDict,
                        "pir" : PirDict,
                        "sff" : SffDict,
                        "swiss" : SwissDict,
                        "qual" : QualDict
                        }

