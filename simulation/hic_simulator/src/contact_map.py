#!/usr/bin/env python
import argparse
import pysam
import numpy as np
import pandas as pd
import sys

# use matplotlib without x-server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

"""
Simulator reads have 'fwd' and 'rev' appended to their pair names
"""
simu_parser = lambda r: (r.qname[:-3], True if r.qname[-3:] == 'fwd' else False)

"""
Generic parse for reads whose paired-ends share the same name. We rely on the status
of the read (read1 or read2).
"""
generic_parser = lambda r: (r.qname, r.is_read1)


class ContactMap:

    def __init__(self, bam, bin_size=50000, simu_reads=False, per_contig=False):
        self.bam = bam
        self.simu_reads = simu_reads
        self.total_seq = len(bam.references)
        self.total_len = sum(bam.lengths)
        self.total_reads = bam.count()
        self.bin_size = bin_size
        self.per_contig = per_contig
        print 'Map based upon mapping containing:\n' \
              '\t{0} sequences\n' \
              '\t{1}bp total length\n' \
              '\t{2} mapped reads'.format(self.total_seq, self.total_len, self.total_reads)
        if per_contig:
            self.bin_count = self.total_seq
            print 'Ignoring bin size since bins are per contig'
            print 'Map details:\n' \
                  '\t{0}bp width\n' \
                  '\t{1}x{1} dimension'.format(self.bin_size, self.bin_count)
        else:
            self.bin_count = int(self.total_len / self.bin_size) + 1
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

    def build_pairs(self):

        # Set the parser behaviour depending on input reads
        if self.simu_reads:
            _parser = simu_parser
        else:
            _parser = generic_parser

        n = 0
        print_rate = self.total_reads if self.total_reads < 1000 else self.total_reads / 1000

        # build a dictionary of all pairs.
        for k, v in self.offsets.iterrows():
            for r in self.bam.fetch(v['name']):
                n += 1
                if n % print_rate == 0:
                    msg = 'Processing {0}/{1} {2}'.format(k+1, self.total_seq, v['name'])
                    progress(n, self.total_reads, msg)

                if r.is_secondary:
                    continue

                rn, rdir = _parser(r)
                if self.per_contig:
                    ix = k
                else:
                    ix = int((v['delta'] + r.pos) / self.bin_size)
                if rn not in self.pairs:
                    self.pairs[rn] = {True: [], False: []}
                self.pairs[rn][rdir].append(ix)

        print '\nFinished building pairs'
        print 'Pairs {0}'.format(len(self.pairs))

    def map_legend(self):
        # TODO
        raise RuntimeError('Unimplemented')

    def _init_map(self, dt=np.int32):
        print 'Initialising contact map of {0}x{0} from total extent of {1}bp over {2} sequences'.format(
            self.bin_count, self.total_len, self.total_seq)
        return np.zeros((self.bin_count, self.bin_count), dtype=dt)

    def calculate_map(self):
        """
        Calculate a raw contact map once .build_pairs() as been called.
        :return: np array representing the contact map
        """
        print 'Beginning calculation of contact map'
        self.raw_map = self._init_map()

        total_pairs = len(self.pairs)
        n = 0
        rate = total_pairs / 1000

        unpaired = 0
        _map = self.raw_map
        for v in self.pairs.values():
            n += 1
            if n % rate == 0:
                progress(n, total_pairs, 'Accumulating')

            if len(v[True]) < 1 or len(v[False]) < 1:
                unpaired += 1
                continue

            for ir in v[True]:
                for ic in v[False]:
                    if ic > ir:
                        _map[ir][ic] += 1
                    else:
                        _map[ic][ir] += 1

        print '\nIgnored {0} unpaired contacts'.format(unpaired)
        print '\nFinished calculation of contact map'
        print 'Total raw map weight {0}'.format(np.sum(_map))

    def calculate_scaled_map(self):
        """
        Calculate a scaled contact map from the raw map. If the raw map has not yet been
        calculated (.calculate_map()) then it will be called first.

        Currently, scaling tries to normalise raw frequencies by contig length and
        read depth.

        :return: np.array representing scaled values.
        """

        if self.norm_map:
            print 'Returning previously calculated normalised map'
            return self.norm_map

        elif self.raw_map is None:
            print 'Raw map has not been calculated and must be calculated first'
            self.calculate_map()

        #
        # TODO need to apply this individually. It will avoid fudge-factor background
        #

        # simple within-contig normalisation
        # this will not adjust weights inter-contigs
        # main intention is to put contigs on similar
        # footing when differing in read-richness and length

        # make a copy of raw contacts, but in double float
        _map = self.raw_map.astype(np.float64)

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

        # global normalisation so map sums to 1.
        _map /= np.sum(_map)
        self.norm_map = _map
        print '\nFinished scaling contact map'

    def calculate_block_map(self):
        """
        Calculate a map where each bin represents an entire contig.
        :return: np.array representing contacts per contig
        """
        cumlen = np.cumsum(self.offsets['length'].as_matrix(), dtype=np.int32)
        end_points = np.column_stack((np.insert(cumlen[:-1], 0, 0), cumlen)) / self.bin_size
        # bins per contig, remove any blocks with zero size
        blocks = end_points[:, 1] - end_points[:, 0]
        blocks = blocks[blocks > 0]
        print blocks

        block_means = np.zeros((blocks.shape[0], blocks.shape[0]))
        n = 0
        for i in xrange(len(blocks)):
            m = n
            for j in xrange(i, len(blocks)):
                print '{0}:{1},{2}:{3}'.format(n, n+blocks[i], m, m+blocks[j])
                # extra block submatrix
                sm = self.raw_map[n:n+blocks[i], m:m+blocks[j]]
                m = np.ma.masked_where(sm == 0, sm).mean()
                if np.isnan(m):
                    print 'Was nan'
                    m = 0.0
                block_means[i, j] = m
                m += blocks[j]
            n += blocks[i]

    def plot_map(self, pname, normalised=False, remove_diag=False):
        """
        Generate a plot (png) for the contact map.
        :param pname: the plot file name
        :param normalised: True - scaled map, False - raw map
        :param remove_diag: True - remove the central diagonal possibly improving the mapping of colour range
        to dynamic range of contact frequencies.
        """
        if normalised:
            # add a minimum value that is twice as small as the smallest non-zero element.
            # this just forms the lower-bound before taking log in plot.
            _map = self.norm_map.copy()
            _map += np.ma.masked_where(_map == 0, _map).min() / 2
        else:
            # being raw counts, many elements are zero
            # a copy is made to avoid side-effect of adding minimum value before taking logs
            _map = self.raw_map.astype(np.float64)
            # being pure counts, we just add 1 to every element. After log transformation
            # the empty elements will be zero.
            _map += 1.0

        if remove_diag:
            np.fill_diagonal(_map, np.min(_map))

        fig = plt.figure(frameon=False)
        img_h = float(self.bin_count) / 100
        fig.set_size_inches(img_h, img_h)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(np.log(_map), interpolation='nearest')
        fig.savefig(pname, dpi=100)

    def write_map(self, oname, normalised=False):
        """
        Write a contact map as an ascii table to a file.
        :param oname: the file name to write
        :param normalised: True - write scaled map, False - write raw map
        """
        _map = self.norm_map if normalised else self.raw_map
        np.savetxt(oname, _map)


def progress(count, total, suffix=''):
    """
    Simple progress indicator for command line.
    """
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create a HiC contact map from mapped reads')
    #parser.add_argument('-q', '--map-qual', default=0, help='Minimum acceptable mapping quality')
    parser.add_argument('--per-contig', default=False, action='store_true', help='Bins are per contig')
    parser.add_argument('--simu-reads', default=False, action='store_true', help='Handle simulator reads')
    parser.add_argument('--bin-size', type=int, default=25000, help='Bin size in bp (25000)')
    parser.add_argument('--remove-diag', default=False, action='store_true', help='Remove the central diagonal from plot')
    parser.add_argument('bamfile', metavar='BAMFILE', nargs=1, help='BAM file to read')
    parser.add_argument('output', metavar='OUTPUT_BASE', nargs=1, help='Output base name')
    args = parser.parse_args()

    with pysam.AlignmentFile(args.bamfile[0], 'rb') as bam:

        contacts = ContactMap(bam, bin_size=args.bin_size, simu_reads=args.simu_reads, per_contig=args.per_contig)
        contacts.build_pairs()
        contacts.calculate_map()

        #contacts.calculate_block_map()

        print 'Writing raw output'
        contacts.write_map('{0}.raw.cm'.format(args.output[0]), normalised=False)
        contacts.plot_map('{0}.raw.png'.format(args.output[0]), normalised=False, remove_diag=args.remove_diag)

        if not args.per_contig:
            contacts.calculate_scaled_map()
            print 'Writing scaled output'
            contacts.write_map('{0}.scl.cm'.format(args.output[0]), normalised=True)
            contacts.plot_map('{0}.scl.png'.format(args.output[0]), normalised=True, remove_diag=args.remove_diag)
