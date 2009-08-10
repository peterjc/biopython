from Bio import SeqIO
#NOTE - Lots of lazy imports further on...

def _genbank_convert(in_handle, in_format, out_handle, out_format, alphabet=None) :
    assert in_format in ["genbank","gb"]
    #We don't need to parser the features...
    from Bio.GenBank.Scanner import GenBankScanner
    records = GenBankScanner().parse_records(in_handle, do_features=False)
    if alphabet :
        records = SeqIO._force_alphabet(records)
    return SeqIO.write(records, out_handle, out_format)

def _embl_convert(in_handle, in_format, out_handle, out_format, alphabet=None) :
    assert in_format == "embl"
    #We don't need to parser the features...
    from Bio.GenBank.Scanner import EmblScanner
    records = EmblScanner().parse_records(in_handle, do_features=False)
    if alphabet :
        records = SeqIO._force_alphabet(records)
    return SeqIO.write(records, out_handle, out_format)


_converter = {
    ("genbank", "fasta") : _genbank_convert,
    ("genbank", "tab") : _genbank_convert,
    ("gb", "fasta") : _genbank_convert,
    ("gb", "tab") : _genbank_convert,
    ("embl", "fasta") : _embl_convert,
    ("embl", "tab") : _embl_convert,
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
    
