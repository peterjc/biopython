#!/usr/bin/env python

import os

h = open("tags.sam", "w")
h.write("@SQ\tSN:chr1\tLN:100\n")
h.write("@SQ\tSN:chr2\tLN:200\n")

#Single integers (SAM code i, BAM codes cCsSiI), sorted by abs
i_values = [0,1,-1,7,8,9,15,16,17,126,127,128,254,255,256,257,-253,-254,-255,-256,
            32000,33000,-33000,64000,-64000,1234567890,-1234567890,
            2147483647, -2147483647]
#Stop here due to bug in samtools 0.1.18 with anything outside uint32
#http://sourceforge.net/mailarchive/message.php?msg_id=28252016
for i in i_values:
    h.write("tag_xx:i:%i\t0\tchr1\t1\t255\t4X\t*\t0\t0\tACGT\t<<<<\txx:i:%i\n" % (i,i))

#Now write some arrays of integers
for code, lower, upper in [("i", -(2**31), 2**31 - 1),
                           ("I",        0, 2**32 - 1)]:
    for i in range(1,len(i_values)+1):
        v = i_values[:i]
        if lower <= min(v) and max(v) <= upper:
            v = ",".join(map(str, i_values[:i]))
            h.write("tag_b%s:B:%s,%s\t0\tchr1\t10\t255\t4X\t*\t0\t0\tACGT\t<<<<\tb%s:B:%s,%s\n" % (code,code,v,code,code,v))
    
#Single precision floats (SAM code f, BAM code f)
f_values = ["0", "-0", "-1.2345", "1.12345", "3.1415e-12", "inf", "-inf", "nan"]
for f in f_values:
    #Allowed [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?
    #These values have been checked by hand to confirm
    #SAM -> BAM -> SAM with samtools preserves them
    #and they match the read name.
    h.write("tag_ff:f:%s\t0\tchr1\t20\t255\t4X\t*\t0\t0\tACGT\t<<<<\tff:f:%s\n" % (f,f))

#Now write it as an array of floats,
for f in [f_values[:5], f_values, f_values[::-1]]:
    v = ",".join(f)
    h.write("tag_bf:B:f,%s\t0\tchr1\t30\t255\t4X\t*\t0\t0\tACGT\t<<<<\tbf:B:f,%s\n" % (v, v))

#Printable character(s) including space (SAM code Z, BAM code Z)
for a in range(32,127):
    #Try both single character (see also code A), and something a little longer.
    for v in [chr(a), chr(a)*3]:
        h.write("tag_zz:Z:%s\t0\tchr1\t40\t255\t4X\t*\t0\t0\tACGT\t<<<<\tzz:Z:%s\n" % (v, v))

#Single printed characters excluding space (SAM code A, BAM code A)
for a in range(32,127):
    v = chr(a)
    h.write("tag_aa:A:%s\t0\tchr1\t50\t255\t4X\t*\t0\t0\tACGT\t<<<<\taa:A:%s\n" % (v, v))

#Hex strings (SAM code H, BAM code H)
for i in [0, 1, 2, 15, 16, 16, 31, 32, 33, 63, 64, 65, 128, 256, 32000, 33000,64000,1234567890]:
    v = "%x" % i
    if len(v) % 2 == 1:
        #Odd length, add leading zero to follow spec
        v = "0" + v
    if v.upper() == v.lower():
        values = [v]
    else:
        #Gets stored as a string in BAM, so can preserve the case
        values = [v.upper(), v.lower()]
    for v in values:
        h.write("tag_hh:H:%s\t0\tchr1\t60\t255\t4X\t*\t0\t0\tACGT\t<<<<\thh:H:%s\n" % (v, v))
h.close()
print "SAM done"

assert 0 == os.system("samtools view -S -b tags.sam > tags.bam")
assert 0 == os.system("samtools index tags.bam")
print "BAM done"
