#!/usr/bin/env python
import argparse
import pysam
import numpy as np
import pandas as pd
import sys

# matplotlib without xserver
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


class ContactMap:
    MIN_CONTACT_VALUE = 0.1

    def __init__(self, bam, bin_size=50000, simu_reads=False):
        self.bam = bam
        self.simu_reads = simu_reads
        self.total_seq = len(bam.references)
        self.total_len = sum(bam.lengths)
        self.total_reads = bam.count()
        self.bin_size = bin_size
        self.bin_count = int(self.total_len / self.bin_size) + 1
        print 'Map based upon mapping containing:\n' \
              '\t{0} sequences\n' \
              '\t{1}bp total length\n' \
              '\t{2} mapped reads'.format(self.total_seq, self.total_len, self.total_reads)
        print 'Map details:\n' \
              '\t{0}bp width\n' \
              '\t{1}x{1} dimension'.format(self.bin_size, self.bin_count)
        self.raw_map = None
        self.norm_map = None

        self.offsets = {}
        for seqname in bam.references:
            tid = bam.gettid(seqname)
            self.offsets[tid] = {'name': seqname,
                                 'delta': sum(bam.lengths[:tid]),
                                 'reads': bam.count(seqname),
                                 'length': bam.lengths[tid]}
        self.offsets = pd.DataFrame.from_dict(self.offsets, orient='index')
        self.offsets.index.name = 'tid'
        self.pairs = {}

    def refname(self, tid):
        return self.bam.references[tid]

    def build_pairs(self):
        if self.simu_reads:
            n = 0
            for k, v in self.offsets.iterrows():
                for r in self.bam.fetch(v['name']):
                    n += 1
                    if n % 100 == 0:
                        msg = 'Processing {0}/{1} {2}'.format(k+1, self.total_seq, v['name'])
                        progress(n, self.total_reads, msg)

                    rn, rdir = r.qname[:-3], r.qname[-3:]
                    ix = int((v['delta'] + r.pos) / self.bin_size)
                    if rn not in self.pairs:
                        self.pairs[rn] = {'fwd': [], 'rev': []}
                    self.pairs[rn][rdir].append(ix)
        else:
            n = 0
            for k, v in self.offsets.iterrows():
                for r in self.bam.fetch(v['name']):
                    n += 1
                    if n % 100 == 0:
                        msg = 'Processing {0}/{1} {2}'.format(k+1, self.total_seq, v['name'])
                        progress(n, self.total_reads, msg)

                    # try to skip useless reads, unpaired/unmapped, etc.
                    # this does not seem to capture all useless
                    if r.is_unmapped or not r.is_paired or r.mate_is_unmapped:
                        continue

                    # simply record all positions, we'll later ignore any
                    # pair with more than 2 locations... dumb!
                    rn = r.qname
                    ix = int((v['delta'] + r.pos) / self.bin_size)
                    if rn not in self.pairs:
                        self.pairs[rn] = [ix]
                    else:
                        self.pairs[rn].append(ix)

        print '\nFinished building pairs'
        print 'Pairs {0}'.format(len(self.pairs))

    def _init_map(self, dt=np.int32):
        print 'Initialising contact map of {0}x{0} from total extent of {1}bp over {2} sequences'.format(
            self.bin_count, self.total_len, self.total_seq)
        return np.zeros((self.bin_count, self.bin_count), dtype=dt)

    def calculate_map(self):
        print 'Beginning calculation of contact map'
        self.raw_map = self._init_map()

        n = 0
        total_pairs = len(self.pairs)
        _map = self.raw_map
        if self.simu_reads:
            for v in self.pairs.values():
                n += 1
                if n % 100 == 0:
                    progress(n, total_pairs, 'Accumulating')

                for ir in v['fwd']:
                    for ic in v['rev']:
                        # we do this to create a upper-tri only matrix
                        if ic > ir:
                            _map[ir][ic] += 1
                        else:
                            _map[ic][ir] += 1
        else:
            skipped = 0
            for v in self.pairs.values():
                n += 1
                if n % 100 == 0:
                    progress(n, total_pairs, 'Accumulating')
                if len(v) == 2:
                    ir, ic = v
                    _map[ir][ic] += 1
                else:
                    skipped += 1
            print '\nIgnored {0} abnormal contacts'.format(skipped)

        print '\nFinished calculation of contact map'
        print 'Total raw map weight {0}'.format(np.sum(_map))

    def calculate_scaled_map(self):

        if self.norm_map:
            print 'Returning previously calculated normalised map'
            return self.norm_map

        elif self.raw_map is None:
            print 'Raw map has not been calculated and must be calculated first'
            self.calculate_map()

        # simple within-contig normalisation
        # this will not adjust weights inter-contigs
        # main intention is to put contigs on similar
        # footing when differing in read-richness and length

        # make a copy of raw contacts, but in double float
        _map = self.raw_map.astype(np.float64)
        _map += ContactMap.MIN_CONTACT_VALUE

        for i in xrange(len(self.offsets)):

            # relative length and number of reads per contig
            rel_reads = self.offsets['reads'][i] / float(self.total_reads)
            rel_length = self.offsets['length'][i] / float(self.total_len)

            # determine submatrix range
            start = int(self.offsets['delta'][i] / self.bin_size)
            if i < len(self.offsets) - 1:
                end = int(self.offsets['delta'][i + 1] / self.bin_size)
            else:
                end = self.bin_count

            # avoid divide by zero
            if rel_reads == 0.0:
                rel_reads = 1.0
            if rel_length == 0.0:
                rel_length = 1.0

            scl = 1.0 / (rel_reads * rel_length)
            #print '{0}/{1} {2} submatrix[{3},{4}][{3},{4}], 1/({5:.3e} * {6:.3e}) = scl {7:.3e}'.format(
            #    i+1, self.total_seq, self.offsets['name'][i], start, end, rel_reads, rel_length, scl)
            progress(i, len(self.offsets), 'Scaling {0} by {1:.2e}'.format(self.offsets['name'][i], scl))

            # apply scale factor only to the upper triangle
            _sm = _map[start:end, start:end]
            _map[start:end, start:end] = np.tril(_sm, -1) + scl * np.triu(_sm)

        _map /= np.sum(_map)
        self.norm_map = _map
        print '\nFinished scaling contact map'

    def plot_map(self, pname, normalised=False):

        if normalised:
            # the version is already prepared for taking log
            _map = self.norm_map
        else:
            # being raw counts, many elements are zero
            # a copy is made to avoid side-effect of adding minimum value before taking logs
            _map = self.raw_map.astype(np.float64)
            _map += ContactMap.MIN_CONTACT_VALUE

        fig = plt.figure(frameon=False)
        img_h = float(self.bin_count) / 100
        fig.set_size_inches(img_h, img_h)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(np.log(_map), interpolation='nearest')
        fig.savefig(pname, dpi=100)

    def write_map(self, oname, normalised=False):
        _map = self.norm_map if normalised else self.raw_map
        np.savetxt(oname, _map)


def progress(count, total, suffix=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))


parser = argparse.ArgumentParser(description='Create a HiC contact map from mapped reads')
parser.add_argument('--simu-reads', default=False, action='store_true', help='Handle simulator reads')
parser.add_argument('--bin-size', type=int, default=25000, help='Bin size in bp (25000)')
parser.add_argument('bamfile', metavar='BAMFILE', nargs=1, help='BAM file to read')
parser.add_argument('output', metavar='OUTPUT_BASE', nargs=1, help='Output base name')
args = parser.parse_args()

with pysam.AlignmentFile(args.bamfile[0], 'rb') as bam:

    contacts = ContactMap(bam, bin_size=args.bin_size, simu_reads=args.simu_reads)
    contacts.build_pairs()
    contacts.calculate_map()
    contacts.calculate_scaled_map()

    # output results
    print 'Writing matrices'
    contacts.write_map('{0}.raw.cm'.format(args.output[0]), normalised=False)
    contacts.write_map('{0}.scl.cm'.format(args.output[0]), normalised=True)
    print 'Writing figures'
    contacts.plot_map('{0}.raw.png'.format(args.output[0]), normalised=False)
    contacts.plot_map('{0}.scl.png'.format(args.output[0]), normalised=True)
