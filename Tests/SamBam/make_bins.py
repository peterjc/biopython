#!/usr/bin/env python

import os
from Bio.Sequencing.SamBam import reg2bin

h = open("bins.sam", "w")

max_len = 2**29 - 1
references = [("tiny", 1000),
              ("small", 2**16),
              ("large", 2**20),
              #("max", max_len),
]

bin_lens = [2**14, 2**17, 2**20, 2**23, 2**226, 2**29]

lengths = [1, 2, 3]
for p in range(1,30):
    x = 2**p
    lengths.extend([x, x+1])
lengths.sort()
print "Trying %i lengths" % len(lengths)

starts = set([1, 2, max_len-1, max_len])
for s in range(1, max_len, 2**13):
    starts.add(s)
assert min(starts) == 1
assert max(starts) == max_len, (max(starts), max_len)
starts = sorted(starts)
print "Trying %i starts" % len(starts)

for ref, ref_len in references:
    h.write("@SQ\tSN:%s\tLN:%i\n" % (ref, ref_len))

for ref, ref_len in references:
    #writing SAM, so start=1 is the very first base (base 0 in BAM)
    bins = set()
    count = 0
    for start in starts:
        for length in lengths:
            end = start + length - 1
            if end > ref_len:
                continue
            bin = reg2bin(start-1, end)
            bins.add(bin)
            if length > 5 or start > 5:
                cigars = ["%iM" % length]
            else:
                cigars = ["%iM" % length, "%iX" % length, "%i=" % length,
                          "1I%iM" % length, "%iM1I" % length, "1D%iM" % length]
                if length > 1:
                    cigars.extend(["1M1I%iM" % (length-1), "1M1D%iM" % (length-1)])
                if length > 2:
                    cigars.extend(["1M1D%iM1I1M" % (length-2), "1M1N%iM" % (length-1)])
            seq = "*"
            flag = 0
            for cigar in cigars:
                count += 1
                name = "%s:r%i:%i..%i:len%i:bin%i:hexbin%s" % \
                       (ref, count, start-1, end, length, bin, hex(bin))
                h.write("%s\t%i\t%s\t%i\t255\t%s\t*\t0\t0\t%s\t*\n" \
                        % (name, flag, ref, start, cigar, seq))
    print "%s len %i, %i reads in %i bins:" % (ref, ref_len, count, len(bins))
    print "(bins %i, ..., %i)" % (min(bins), max(bins))
h.close()
print "SAM done"

assert 0 == os.system("samtools view -S -b bins.sam > bins.bam")
assert 0 == os.system("samtools index bins.bam")
print "BAM done"

