import pysam

samfile = pysam.Samfile("merged.sorted.bam", "rb" )
for pileupcolumn in samfile.pileup( 'K03455|HIVHXB2CG', 138, 139):
    print
    print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
    for pileupread in pileupcolumn.pileups:
        print '\tbase in read %s = %s' % (pileupread.alignment.qname, pileupread.alignment.seq[pileupread.qpos])

samfile.close()
