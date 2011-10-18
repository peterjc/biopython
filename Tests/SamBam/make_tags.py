import os

h = open("tags.sam", "w")
h.write("@SQ\tSN:chr1\tLN:100\n")
h.write("@SQ\tSN:chr2\tLN:200\n")
for i in [0,1,-1,126,127,128,254,255,256,257,-253,-254,-255,-256,32000,33000,-33000,64000,64000]:
    h.write("tag_xx:i:%i\t0\tchr1\t50\t255\t4X\t*\t0\t0\tACGT\t<<<<\txx:i:%i\n" % (i,i))
h.close()
print "SAM done"

assert 0 == os.system("samtools view -S -b tags.sam > tags.bam")
print "BAM done"
