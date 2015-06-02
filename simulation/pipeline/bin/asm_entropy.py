#!/usr/bin/env python
# coding=utf-8

import pysam
import numpy as np
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='Calculate assembly weighted entropy from reads to contig BAM')
parser.add_argument('bam_file', metavar='BAM', nargs=1, help='BAM file of Reads to Contigs')
parser.add_argument('--plot', action='store_true', default=False, help='Plot histogram of contig entropy')
args = parser.parse_args()

# determine the source purity of each contig
with pysam.AlignmentFile(args.bam_file[0], 'rb') as sam_file:
    extent = np.sum(sam_file.lengths, dtype=np.float64)
    weight = dict([(rn, sam_file.lengths[sam_file.gettid(rn)] / extent) for rn in sam_file.references])

    purity = dict.fromkeys(sam_file.references)
    for rd in sam_file.fetch():
        ctg_name = sam_file.getrname(rd.reference_id)
        if purity[ctg_name] is None:
            # assume minimum of one read, this avoid zeros in calculation
            purity[ctg_name] = {'A': 1, 'B': 1, 'C': 1, 'D': 1}
        purity[ctg_name][rd.query_name[0]] += 1

# calculate the entropy for each contig, where a contig of
# highly pure source reads will have minimum S.
ctg_entropy = {}
for ctg, rd_counts in purity.iteritems():
    if rd_counts is None:
        print '{0} had no read information'.format(ctg)
    else:
        tot_rd = float(sum(rd_counts.values()))
        xi = np.array(rd_counts.values())
        pi = xi/tot_rd
        ent = -np.sum(pi * np.log2(pi))
        ctg_entropy[ctg] = ent

asm_entropy = 0
for ctg, Si in ctg_entropy.iteritems():
    asm_entropy += weight[ctg] * Si

print 'Weighted assembly entropy: {0}'.format(asm_entropy)

# plot
if args.plot:
    p = plt.hist(ctg_entropy.values(), bins=100)
    plt.axvline(x=asm_entropy, color='r')
    plt.show()
