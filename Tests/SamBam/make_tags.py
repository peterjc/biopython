import os

h = open("tags.sam", "w")
h.write("@SQ\tSN:chr1\tLN:100\n")
h.write("@SQ\tSN:chr2\tLN:200\n")
#Single integers (SAM code i, BAM codes cCsSiI
for i in [0,1,-1,126,127,128,254,255,256,257,-253,-254,-255,-256,
          32000,33000,-33000,64000,-64000,1234567890,-1234567890]:
    h.write("tag_xx:i:%i\t0\tchr1\t50\t255\t4X\t*\t0\t0\tACGT\t<<<<\txx:i:%i\n" % (i,i))
#Printable character(s) including space (SAM code Z, BAM code Z)
for a in range(32,127):
    #Try both single character (see also code A), and something a little longer.
    for v in [chr(a), chr(a)*3]:
        h.write("tag_zz:Z:%s\t0\tchr1\t50\t255\t4X\t*\t0\t0\tACGT\t<<<<\tzz:Z:%s\n" % (v, v))
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
        h.write("tag_hh:H:%s\t0\tchr1\t50\t255\t4X\t*\t0\t0\tACGT\t<<<<\thh:H:%s\n" % (v, v))
h.close()
print "SAM done"

assert 0 == os.system("samtools view -S -b tags.sam > tags.bam")
print "BAM done"
