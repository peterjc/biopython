#!/usr/bin/env python
# Copyright 2012 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
usage = """Python script using Biopython to mimic 'samtools idxstats'

This script uses the Bio.Sequencing.SamBam module of Biopython to
reimplement the 'samtools idxstats' command in pure Python. This
may prove useful on platforms where samtools is not available,
but its primarily usecase was to allow testing of the BAI file
(BAM index) parsing in Biopython.

Usage:

$ sambam_idxstats.py <in.bam>

i.e. Takes a single command line parameter, the name of a BAM
file. It will look for a matching BAI file under the name of
<in.bam.bai>

Output:

Prints a four column tab separated table to stdout. There is one
line per reference sequence in the BAM file giving the reference
sequence name, reference sequence length (base pairs), number of
reads mapped to that reference, and number of unmapped reads
assigned to that reference (usually from pair-end reads). Finally
the last line gives the unmapped and unplaced reads (with dummy
entries of asterisk, zero and zero for the first three columns).
"""

import sys
import os
import gzip


from Bio.Sequencing import SamBam
from Bio.Sequencing.SamBam import bai

if len(sys.argv) != 2:
    sts.stderr.write("ERROR - Missing required argument (BAM filename)\n\n")
    sys.stderr.write(usage)
    sys.exit(1)
bam_filename = sys.argv[1]

if not os.path.isfile(bam_filename):
    sys.stderr.write("File %s not found\n" % bam_filename)
    sys.exit(1)

bai_filename = bam_filename + ".bai"
if not os.path.isfile(bai_filename):
    sys.stderr.write("BAM index file %s not found\n" % bai_filename)
    sys.exit(1)
#TODO - Try replacing .bam with .bai instead

def idxstats(bam_filename, bai_filename):
    #Don't need random access, so can just use gzip not bgzf
    handle = gzip.open(bam_filename, "rb")
    header_text, ref_count = SamBam._bam_file_header(handle)
    references = [SamBam._bam_file_reference(handle) for i in range(ref_count)]
    handle.close()

    handle = open(bai_filename, "rb")
    indexes, unmapped = bai._load_bai(handle)
    if unmapped is None:
        sys.stderr.write("ERROR - Old index lacks unmapped read information, re-index your BAM file\n")
        sys.exit(1)

    if len(indexes) != len(references):
        sys.write("ERROR: BAM file has %i references, BAI has %i\n" \
                  % (len(references), len(indexes)))
        sys.exit(1)

    for (reference, length), (chunks, linear, mapped, ref_unmapped, u_start, u_end) in zip(references, indexes):
        if mapped is None:
            mapped = 0
        if ref_unmapped is None:
            ref_unmapped = 0
        #TODO - Check length versus linear bins
        yield reference, length, mapped, ref_unmapped
    yield "*", 0, 0, unmapped

for name, length, mapped, placed in idxstats(bam_filename, bai_filename):
    print "%s\t%s\t%i\t%i" % (name, length, mapped, placed)
