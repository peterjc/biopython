# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Optimised sequence conversion code (PRIVATE).

You are not expected to access this module, or any of its code, directly. This
is all handled internally by the Bio.SeqIO.convert(...) function which is the
public interface for this.

The idea here is rather while doing this will work:

from Bio import SeqIO
records = SeqIO.parse(in_handle, in_format)
count = SeqIO.write(records, out_handle, out_format)

it is shorter to write:

from Bio import SeqIO
count = SeqIO.convert(in_handle, in_format, out_handle, out_format)

Also, the convert function can take a number of special case optimisations. This
means that using Bio.SeqIO.convert() may be faster, as well as more convenient.
"""

from Bio import SeqIO
#NOTE - Lots of lazy imports further on...

def _genbank_convert(in_handle, in_format, out_handle, out_format, alphabet=None) :
    """Fast GenBank to "other" conversion, where other does not use features (PRIVATE)."""
    assert in_format in ["genbank","gb"]
    #We don't need to parser the features...
    from Bio.GenBank.Scanner import GenBankScanner
    records = GenBankScanner().parse_records(in_handle, do_features=False)
    if alphabet :
        records = SeqIO._force_alphabet(records)
    return SeqIO.write(records, out_handle, out_format)

def _embl_convert(in_handle, in_format, out_handle, out_format, alphabet=None) :
    """Fast EMBL to "other" conversion, where other does not use features (PRIVATE)."""
    assert in_format == "embl"
    #We don't need to parser the features...
    from Bio.GenBank.Scanner import EmblScanner
    records = EmblScanner().parse_records(in_handle, do_features=False)
    if alphabet :
        records = SeqIO._force_alphabet(records)
    return SeqIO.write(records, out_handle, out_format)

def _fastq_convert_fasta(in_handle, in_format, out_handle, out_format, alphabet=None) :
    """Fast FASTQ to FASTA conversion (PRIVATE).

    Avoids dealing with the FASTQ quality encoding, and creating SeqRecord and
    Seq objects in order to speed up this conversion. For example,

    >>> expected = open("Quality/example.fasta", "ru").read()
    >>> from StringIO import StringIO
    >>> handle = StringIO("")
    >>> converted = _fastq_convert_fasta(open("Quality/example.fastq", "rU"),
    ...                                  "fastq", handle, "fasta")
    >>> expected == handle.getvalue()
    True
    """
    assert in_format in ["fastq", "fastq-sanger","fastq-solexa", "fastq-illumina"]
    assert out_format == "fasta"
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    #For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    for title, seq, qual in FastqGeneralIterator(in_handle) :
        count += 1
        out_handle.write(">%s\n" % title)
        #Do line wrapping
        for i in range(0, len(seq), 60):
            out_handle.write(seq[i:i+60] + "\n")
    return count

def _fastq_convert_tab(in_handle, in_format, out_handle, out_format, alphabet=None) :
    """Fast FASTQ to simple tabbed conversion (PRIVATE)."""
    assert in_format in ["fastq", "fastq-sanger","fastq-solexa", "fastq-illumina"]
    assert out_format == "tab"
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    #For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    for title, seq, qual in FastqGeneralIterator(in_handle) :
        count += 1
        out_handle.write("%s\t%s\n" % (title, seq))
    return count

_converter = {
    ("genbank", "fasta") : _genbank_convert,
    ("genbank", "tab") : _genbank_convert,
    ("gb", "fasta") : _genbank_convert,
    ("gb", "tab") : _genbank_convert,
    ("embl", "fasta") : _embl_convert,
    ("embl", "tab") : _embl_convert,
    ("fastq", "fasta") : _fastq_convert_fasta,
    ("fastq-sanger", "fasta") : _fastq_convert_fasta,
    ("fastq-solexa", "fasta") : _fastq_convert_fasta,
    ("fastq-illumina", "fasta") : _fastq_convert_fasta,
    ("fastq", "tab") : _fastq_convert_tab,
    ("fastq-sanger", "tab") : _fastq_convert_tab,
    ("fastq-solexa", "tab") : _fastq_convert_tab,
    ("fastq-illumina", "tab") : _fastq_convert_tab,
    }

def _handle_convert(in_handle, in_format, out_handle, out_format, alphabet=None):
    """SeqIO conversion function (PRIVATE)."""
    try :
        f = _converter[(in_format, out_format)]
    except KeyError :
        f = None
    if f :
        return f(in_handle, in_format, out_handle, out_format, alphabet)
    else :
        records = SeqIO.parse(in_handle, in_format, alphabet)
        return SeqIO.write(records, out_handle, out_format)
    

def _test():
    """Run the Bio.SeqIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os
    if os.path.isdir(os.path.join("..","..","Tests")) :
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join("..","..","Tests"))
        assert os.path.isfile("Quality/example.fastq")
        assert os.path.isfile("Quality/example.fasta")
        doctest.testmod(verbose=0)
        os.chdir(cur_dir)
        del cur_dir
        print "Done"
        
if __name__ == "__main__" :
    _test()
