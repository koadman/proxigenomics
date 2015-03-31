#!/usr/bin/env python
import sys
import pysam
import argparse

parser = argparse.ArgumentParser(description='Filter BAM/SAM files for read coverage and mapping quality')
parser.add_argument('-s', dest='sam_input', default=False, action='store_true', help='Input is ASCII SAM format')
parser.add_argument('--mincov', type=float, default=0.95, help='Minimum coverage of aligned read [0..1]')
parser.add_argument('--minqual', type=int, default=40, help='Minimum mapping quality score of aligned read [40]')
parser.add_argument('input', metavar='INPUT_FILE', nargs=1, help='Input BAM/SAM file')
parser.add_argument('output', metavar='OUTPUT_FILE', nargs=1, help='Output BAM file')
args = parser.parse_args()

# set these here to eliminate some deferencing in the callback
min_cov = args.mincov
min_mapq = args.minqual

# Set input mode based commandline option
input_mode = 'r' if args.sam_input else 'rb'
infile = pysam.AlignmentFile(args.input[0], input_mode)
outfile = pysam.AlignmentFile(args.output[0], 'wb', template=infile)

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
