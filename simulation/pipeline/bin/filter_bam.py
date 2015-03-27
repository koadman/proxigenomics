#!/usr/bin/env python
import sys
import pysam


if len(sys.argv) != 5:
    print 'Usage: [mcov] [mapq] [in file] [out bam]'
    sys.exit(1)

min_cov = float(sys.argv[1])
min_mapq = int(sys.argv[2])

infile = pysam.AlignmentFile(sys.argv[3], 'rb')
outfile = pysam.AlignmentFile(sys.argv[4], 'wb', template=infile)

kept = 0
total = 0


def keepers_callback(read):
    """
    Exclude reads which do not meet the minimum coverage and mapping quality.
    Write successful reads to output bam file.

    Coverage is determined by directly counting up matching regions in the
    CIGAR string and comparing that to the total extent reported in the
    CIGAR.

    Note: this should be a callback, but callbacks do not appear to be invokved
    with fetch in pysam -- despite the documentation saying otherwise.

    :param read: pysam.AlignedSegment
    """
    global kept, total, outfile
    if read.cigarstring is not None:
        mlen = float(sum([cig_i[1] for cig_i in read.cigartuples if cig_i[0] == 0]))
        rlen = float(sum([cig_i[1] for cig_i in read.cigartuples]))
        if mlen/rlen >= min_cov and read.mapping_quality >= min_mapq:
            kept += 1
            outfile.write(read)
    total += 1

for read in infile.fetch():
    keepers_callback(read)

infile.close()
outfile.close()

print 'Kept {0} of {1} {2:.1f}%'.format(kept, total, 100.0*kept/total)
