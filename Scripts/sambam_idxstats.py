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

from Bio.Sequencing.SamBam.bai import idxstats

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

for name, length, mapped, placed in idxstats(bam_filename, bai_filename):
    print "%s\t%s\t%i\t%i" % (name, length, mapped, placed)

