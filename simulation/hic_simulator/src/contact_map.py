#!/usr/bin/env python
import pysam
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def get_pairs(sname):
    """
    Build a dictionary keyed by insert name of mapped positions
    :param sname: sequence name
    :return: dictionary of insert positions
    """
    pairs = {}
    for aln in bam_file.fetch(sname):
        nm, dir = aln.qname[:-3], aln.qname[-3:]
        if nm not in pairs:
            pairs[nm] = {}
        pairs[nm][dir] = aln.pos
    return pairs

parser = argparse.ArgumentParser(description='Create a HiC contact map from mapped reads')
parser.add_argument('--plot', default=False, action='store_true', help='Plot resulting matrix')
parser.add_argument('--bin-size', type=int, default=25000, help='Bin size in bp (25000)')
parser.add_argument('--pseudo-bidir', default=False, action='store_true', help='Produce pseudo-bidirectional distances')
parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
parser.add_argument('bamfile', metavar='BAMFILE', nargs=1, help='BAM file to read')
parser.add_argument('output', metavar='OUTPUT', type=argparse.FileType('w'), help='Output (stdout)')
args = parser.parse_args()

with pysam.AlignmentFile(args.bamfile[0]) as bam_file:

    bin_size = args.bin_size
    total_seq = len(bam_file.lengths)
    total_len = sum(bam_file.lengths)
    bin_count = int(total_len/bin_size) + 1
    contacts = np.zeros((bin_count, bin_count), dtype=np.int32)
    first_bin = 0

    for n, sname in enumerate(bam_file.references):

        if args.verbose and n % 50 == 0:
            print 'Processing {0}/{1} {2:.1f}% - {3}'.format(n, total_seq, n*100.0/total_seq, sname)

        pairs = get_pairs(sname)
        ix = ['fwd','rev']
        for pk, vals in pairs.iteritems():
            if len(vals) != 2:
                continue

            # this just creates the pseudo-bidirectional reads we don't make with simForward.py
            if args.pseudo_bidir:
                np.random.shuffle(ix)

            ir = int((vals[ix[0]] + first_bin) / bin_size)
            ic = int((vals[ix[1]] + first_bin) / bin_size)
            contacts[ir][ic] += 1

        first_bin += bam_file.lengths[n]

    np.savetxt(args.output, contacts)

    if args.plot:
        im = plt.imshow(np.log(contacts), interpolation='nearest')
        plt.show()
