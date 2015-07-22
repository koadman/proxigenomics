#!/usr/bin/env python
import argparse
import pysam
import numpy as np
import math
import sys
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch
from scipy.misc import factorial
from scipy.stats import poisson
from scipy.stats import geom

# use matplotlib without x-server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


import yaml
from collections import OrderedDict
def order_rep(dumper, data):
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', data.items(), flow_style=False)
yaml.add_representer(OrderedDict, order_rep)


def dump_pairs(filename, pairs):
    print 'Dumping pairs to {0}'.format(filename)
    with open(filename, 'w') as out_h:
        # convert numpy int32 to int, yaml chokes otherwise
        ndict = {}
        for rn, reads in pairs.iteritems():
            ndict[rn] = {}
            for k, places in reads.iteritems():
                ndict[rn][k] = map(int, places)  # non-py3
        od = OrderedDict(sorted(ndict.items()))
        yaml.dump(od, out_h)


"""
Simulator reads have 'fwd' and 'rev' appended to their pair names
"""
simu_parser = lambda r: (r.qname[:-3], True if r.qname[-3:] == 'fwd' else False)

"""
Generic parse for reads whose paired-ends share the same name. We rely on the status
of the read (read1 or read2).
"""
generic_parser = lambda r: (r.qname, r.is_read1)


class CutSite(int):
    """
    A cut site is just an annotated integer.
    """
    def __new__(cls, name=None, pos=0):
        i = int.__new__(cls, pos)
        i.name = name
        return i

    def __repr__(self):
        return repr((self.name, int.__repr__(self)))

    def is_same_enzyme(self, other):
        """
        Check two sites share the same enzyme.

        :param other: other site to compare
        :return: True - they are the same (enzyme name)
        """
        return self.name == other.name

    def gap(self, other):
        """
        Absolute separate distance in base-pairs
        :param other: to compare
        :return: absolute value in base-pairs
        """
        return abs(self - other + 1)


class DigestedSequence:

    def __init__(self, enzymes, sequence, is_linear=True):
        self.enzymes = enzymes
        self.sequence = sequence
        self.res_batch = RestrictionBatch(enzymes)
        self.is_linear = is_linear
        self.site_dict = self.res_batch.search(self.sequence.seq, is_linear)

    def get_sites(self):
        """
        Return the set of sites for a given contig, ordered by increasing position.
        :return: list of CutSites
        """
        cutSites = []
        for e_name, ctg_locs in self.site_dict.iteritems():
            for loc in ctg_locs:
                cutSites.append(CutSite(e_name, loc))

        return sorted(cutSites)

    def get_fragments(self):
        """
        Return the genomic fragments resulting from the digestion.
        :return: list of SeqRecords.
        """
        sites = self.get_sites()
        seq = self.sequence

        if len(sites) == 0:
            return seq

        frags = []
        for idx in xrange(1, len(sites)):
            a = sites[idx - 1]
            b = sites[idx]
            frg = seq[a:b]
            frg.id = '{0}:{1}:{2}'.format(frg.id, a, b)
            frg.name = frg.id
            frg.description = 'restriction digest fragment from {0} to {1}'.format(a, b)
            frags.append(frg)

        return frags

    @staticmethod
    def digestion_sites(seq_list, enzyme_names=[]):
        """
        Return a list of sites per sequence, preserving the input list order.
        :param seq_list: list of sequences to analyze
        :param enzyme_names: enzyme used in digestion
        :return: list of sites per sequence
        """
        sites = []
        for seq in seq_list:
            ds = DigestedSequence(enzyme_names, seq)
            seq_sites = ds.get_sites()
            if len(seq_sites) == 0:
                print '\tSequence {0} (length {1}) contained zero sites'.format(seq.id, len(seq))
            sites.append({'name': seq.id, 'pos': seq_sites})
        return sites


class Grouping:

    # masked bit
    MASK = -1

    def __init__(self, site_list, bin_width):
        total_sites = 0
        self.bins = []
        self.map = []
        for seq_sites in site_list:

            if len(seq_sites['pos']) == 0:
                # Sequences which had no sites get empty place holders
                # this maintains order of lists
                self.bins.append(Grouping.MASK)
                self.map.append(None)
                continue

            max_n = len(seq_sites['pos']) - 1
            total_sites += max_n
            bin_n = int(math.floor(max_n / bin_width))

            if bin_n == 0:
                # sequences with <= bin_width sites get one bin
                print '\tSequence {0} had only {1} sites'.format(seq_sites['name'], len(seq_sites['pos']))
                bin_n = 1

            print '\tHighest index {0} requires {1} bins, scale={2:.4f}'.format(max_n, bin_n, bin_n/float(max_n))
            scaled_indices = np.array(range(max_n+1)) * bin_n/float(max_n)
            grp = np.digitize(scaled_indices, np.arange(bin_n)) - 1
            self.bins.append(bin_n)
            self.map.append(np.column_stack((seq_sites['pos'], grp)))
            print '\tGroup sizes {0}'.format(np.bincount(grp))

        self.total = total_sites

        # calculate centers of each fragment groupings. This will be used
        # when measuring separation.
        self.centers = []
        for n, num_bins in enumerate(self.bins):
            if num_bins == -1:
                self.centers.append(None)
            else:
                self.centers.append(np.array([
                    np.mean(self.map[n][np.where(self.map[n][:, 1] == i)]) for i in xrange(num_bins)]))

    def total_bins(self):
        mask_bins = np.ma.masked_equal(self.bins, Grouping.MASK)
        return np.ma.sum(mask_bins)

    def find_nearest(self, ctg_idx, x, above=True):
        # TODO this is very slow! Masked arrays have a severe penalty, therefore
        # need to do this another way. The bounds checking here is no faster than
        # testing the masked array for "all masked".
        group_map = self.map[ctg_idx]
        if above:
            if x > group_map[-1, 0]:
                raise RuntimeError('No site above {0}, largest was {1}'.format(x, group_map[-1, 0]))
            diff = group_map[:, 0] - x
        else:
            if x < group_map[0, 0]:
                raise RuntimeError('No site before {0}, smallest was {1}'.format(x, group_map[0, 0]))
            diff = x - group_map[:, 0]

        # mask any value less than zero and then find the smallest argument.
        masked_diff = np.ma.masked_less(diff, 0)
        idx = masked_diff.argmin()
        return group_map[idx, :]


class SeqOrder:

    def __init__(self, seqs, sites):
        """
        Order is initially determined by the order of supplied sequence list. Indices of sequences
        which possessed no restriction sites are indicated in "no_sites".

        :param seqs: list of sequences
        :param sites: list of re-sites per sequence.
        """
        self.names = [s.id for s in seqs]
        self.lengths = np.array([len(s) for s in seqs])
        self.order = np.array(range(len(seqs)))
        self.no_sites = np.array([len(s['pos']) == 0 for s in sites])

    def is_first(self, a, b):
        """
        Test if a comes before another sequence in order
        :param a: target sequence in order
        :param b: another sequence in order
        :return: True if a comes before b
        """
        ia = np.argwhere(self.order == a)
        ib = np.argwhere(self.order == b)
        return ia < ib

    def intervening(self, a, b):
        """
        For the current internal order, calculate the length of intervening
        sequence between sequence A and sequence B.

        :param a: first sequence
        :param b: second sequence
        :return: total length of sequences currently between A and B in the order.
        """
        ia = np.argwhere(self.order == a)
        ib = np.argwhere(self.order == b)
        if ia > ib:
            ia, ib = ib, ia
        return np.sum(self.lengths[self.order[ia+1:ib]])

    @staticmethod
    def initial_ordering_from_file(fasta_file):
        with open(fasta_file, 'r') as fasta_h:
            seqs = [seq for seq in SeqIO.parse(fasta_h, 'fasta')]
            return SeqOrder(seqs)


class FragmentMap:

    def __init__(self, bam_file, fasta_file, enzyme, bin_width=3, simu_reads=False):

        with open(fasta_file, 'r') as fasta_h:
            seq_list = [seq for seq in SeqIO.parse(fasta_h, 'fasta')]

            print 'Digesting supplied sequences with {0} ...'.format(enzyme)
            self.sites = DigestedSequence.digestion_sites(seq_list, [enzyme])

            print 'Initializing sequence order ...'
            self.order = SeqOrder(seq_list, self.sites)

        with pysam.AlignmentFile(bam_file, 'rb') as self.bam:
            self.simu_reads = simu_reads
            self.total_seq = len(self.bam.references)
            self.total_len = sum(self.bam.lengths)
            self.total_reads = self.bam.count()
            self.bin_width = bin_width

            print 'Map based upon mapping containing:\n' \
                  '\t{0} sequences\n' \
                  '\t{1}bp total length\n' \
                  '\t{2} mapped reads'.format(self.total_seq, self.total_len, self.total_reads)

            print 'Binning restriction sites into groups of approximately {0} ...'.format(bin_width)
            self.groupings = Grouping(self.sites, self.bin_width)

            print 'Map details:\n' \
                  '\t{0} adjacent fragments per bin\n' \
                  '\t{1} total sites\n' \
                  '\t{2}x{2} matrix'.format(self.bin_width, self.groupings.total_bins(), self.groupings.total_bins())

            self.raw_map = None
            self.norm_map = None
            self.pairs = {}

            # process open BAM file and initialise pairs dictionary
            self._build_pairs()

            # using initialized pairs, populate the contact matrix
            self._calculate_map()

    def set_order(self, new_order):
        self.order.order = np.array(new_order)

    def shuffle_order(self):
        """
        Randomly shuffle the internal order of sequences, but do not include those sequences
        which have no sites.
        """
        _order = self.order.order
        _bins = np.array(self.groupings.bins)

        # select those sequences which have bins
        shuf = _order[_bins != Grouping.MASK]

        # shuffle the subset, which is done in place.
        np.random.shuffle(shuf)

        # return the now shuffled elements to the original
        # array, skipping over those which were masked out.
        n = 0
        num_seq = _order.shape[0]
        for i in xrange(num_seq):
            if _bins[i] == Grouping.MASK:
                continue
            _order[i] = shuf[n]
            n += 1

    def get_submatrix(self, i, j):
        """
        Return the submatrix of the raw contact counts pertaining to the interaction
        between sequence i and j.
        :param i: first sequence id
        :param j: second sequence id
        :return: submatrix (nparray) containing raw counts between these sequences
        """
        _bins = self.groupings.bins
        oi = sum(_bins[0:i])
        oj = sum(_bins[0:j])
        return self.raw_map[oi:oi+_bins[i], oj:oj+_bins[j]]

    def calc_likelihood(self, p0, s0, b):
        """
        Calculate the logLikelihood of a given sequence configuration. The model is adapted from
        GRAAL. Counts are Poisson, with lambda parameter dependent on expected observed contact
        rate as a function of inter-fragment separation. This has been shown experimentally to be
        modelled effectively by a power-law (used here).

        :param p0: minimum probability of a contact
        :param s0: characteristic distance at which probability = p0
        :param b: exponential parameter
        :return:
        """
        Nd = np.sum(self.raw_map)
        sumL = 0.0

        for i in xrange(self.total_seq):

            if self.order.no_sites[i]:
                continue

            for j in xrange(i+1, self.total_seq):

                if self.order.no_sites[j]:
                    continue

                L = self.order.intervening(i, j)
                #print '{0},{1} intervening distance {2}'.format(i, j, L)

                centers_i = self.groupings.centers[i]
                centers_j = self.groupings.centers[j]

                # determine relative origin for measuring separation
                # between sequences. If i comes before j, then distances
                # to j will be measured from the end of i -- and visa versa
                if self.order.is_first(i, j):
                    #print 'i={0} comes before j={1}, measuring from end of i to start of j'.format(i,j)
                    s_i = self.order.lengths[i] - centers_i
                    s_j = centers_j
                else:
                    #print 'i={0} comes after j={1}, measuring from end of j to start of i'.format(i,j)
                    s_i = centers_i
                    s_j = self.order.lengths[j] - centers_j

                # go crazy
                d_ij = np.abs(L + s_i[:, np.newaxis] - s_j)
                #print 'd_ij dimensions {0}'.format(d_ij.shape)
                #print d_ij

                # Here we are using a peice-wise continuous function defined in GRAAL
                # to relate separation distance to Poisson rate parameter mu.
                # Above a certain separation distance, the probability reaches a constant minimum value.
                #q_ij = np.piecewise(d_ij, [d_ij <= s0, d_ij > s0], [lambda s: p0 * (s / s0)**b, lambda s: p0])
                #q_ij = np.piecewise(d_ij, [d_ij <= s0, d_ij > s0], [lambda s: geom.pmf(s, p=6.0e-6), lambda s: 2.3e-6])

                # q_ij = np.piecewise(d_ij, [d_ij<=600000., (d_ij>600000) & (d_ij<=3e6), d_ij>3e6],
                #                     [lambda x: 4.0e-6 * (1. - 4.0e-6)**(x - 1.),
                #                      lambda x: 3.0e-6,
                #                      lambda x: 1.0e-8])
                q_ij = np.piecewise(d_ij, [d_ij<3e6, d_ij>3e6],
                                    [lambda x: 0.5*(1.0 / 3e6 - (1. - 6e-6)**x * np.log(1.0 - 6e-6)),
                                     lambda x: 1.0e-8])

                #print q_ij
                #print 'q_ij dimensions {0}'.format(q_ij.shape)
                n_ij = self.get_submatrix(i, j)
                #print 'n_ij dimensions {0}'.format(n_ij.shape)
                #print n_ij
                # TODO just calculate log likelihood here rather than take log after?
                #p_ij = (Nd * q_ij) ** n_ij / factorial(n_ij) * np.exp(-Nd * q_ij)

                p_ij = poisson.logpmf(n_ij, mu=(Nd * q_ij))
                #print p_ij

                #li = np.sum(np.log(p_ij))
                li = np.sum(p_ij)
                #print '{0},{1} logL {2}'.format(i, j, li)
                #print '{0},{1}'.format(self.order.names[i], self.order.names[j])
                sumL += li
        return sumL

    def _build_pairs(self):

        print 'Building pairs ...'

        # Set the parser behaviour depending on input reads
        if self.simu_reads:
            _parser = simu_parser
        else:
            _parser = generic_parser

        skipped = 0
        n = 0
        print_rate = self.total_reads if self.total_reads < 1000 else self.total_reads / 1000
        bin_offset = 0
        # build a dictionary of all pairs.
        for ctg_idx, seq_name in enumerate(self.order.names):

            if self.groupings.bins[ctg_idx] <= 0:
                # don't bother trying to bin reads for sequences with no sites
                continue

            for r in self.bam.fetch(seq_name):
                n += 1
                if n % print_rate == 0:
                    msg = 'Processing {0}/{1} {2}'.format(ctg_idx+1, self.total_seq, seq_name)
                    progress(n, self.total_reads, msg)

                if r.is_secondary:
                    continue

                rn, rdir = _parser(r)

                # beginning of read is depends in mapping direction
                refpos = r.reference_start if rdir else r.reference_end

                try:
                    bin_i = self.groupings.find_nearest(ctg_idx, refpos, rdir)
                    if rn not in self.pairs:
                        self.pairs[rn] = {True: [], False: []}
                    self.pairs[rn][rdir].append(bin_i[1] + bin_offset)
                except RuntimeError as er:
                    skipped += 1
                    #print '\tSkipping read {0} mapped to {1}: {2}'.format(rn, seq_name, er)

            bin_offset += self.groupings.bins[ctg_idx]
            #print 'ctg_idx {0} bin_offset {1}'.format(ctg_idx, bin_offset)

        print '\nFound {0} fragments in BAM, {1} not reconciled with any site'.format(len(self.pairs), skipped)
        print '\nFinished building pairs'

    def _init_map(self, dt=np.int32):
        n_bins = self.groupings.total_bins()
        print '\tInitialising contact map of {0}x{0} fragment bins, ' \
              'representing {1} bp over {2} sequences'.format(n_bins, self.total_len, self.total_seq)
        return np.zeros((n_bins, n_bins), dtype=dt)

    def _determine_block_shifts(self):
        """
        For the present ordering, calculate the block (whole-group) shifts necessary
        to permute the contact matrix to match the order.
        :return: list of tuples (start, stop, shift)
        """
        # make sure the two arrays are np arrays for fancy indexing tricks
        _order = np.array(self.order.order)
        _bins = np.array(self.groupings.bins)
        _shuf_bins = _bins[_order]

        shifts = []
        # iterate through current order, determine necessary block shifts
        curr_bin = 0
        for shuff_i, org_i in enumerate(_order):
            if _bins[org_i] == Grouping.MASK:
                # skip masked bins
                continue
            intervening_bins = _bins[:org_i]
            shft = -(curr_bin - np.sum(intervening_bins))
            start = np.sum(_bins[:org_i]) if org_i > 0 else 0
            #stop = start+_bins[org_i]
            #print '{0} -> {1} {2} {3},{4} s:{5} {6},{7}'.format(org_i, shuff_i, intervening_bins, start, stop, shft, curr_bin, curr_bin+_bins[org_i])
            shifts.append((curr_bin, curr_bin+_bins[org_i], shft))
            curr_bin += _shuf_bins[shuff_i]
        return shifts

    def _make_permutation_matrix(self):
        """
        Create the permutation matrix required to reorder the contact matrix to
        match the present ordering.
        :return: permutation matrix
        """
        # as a permutation matrix, the identity matrix causes no change
        P = np.identity(self.groupings.total_bins())
        block_shifts = self._determine_block_shifts()
        for si in block_shifts:
            if si[2] == 0:
                # ignore blocks which are zero shift.
                # as we started with identity, these are
                # already defined.
                continue
            # roll along columns, those rows matching each block with
            # a shift.
            pi = np.roll(P, si[2], 1)[si[0]:si[1], :]
            # put the result back into P, overwriting previous identity diagonal
            P[si[0]:si[1], :] = pi
        return P

    def reorder_map(self):
        """
        Reorder the contact matrix to reflect the present order.
        :return: reordered contact matrix
        """
        # this may not be necessary, but we construct the full symmetrical map
        # before reordering, just in case something is reflected over the diagonal
        # and ends up copying empty space
        full_map = np.tril(self.raw_map.transpose(), -1) + self.raw_map
        P = self._make_permutation_matrix()
        # two applications of P is required to fully reorder the matrix
        # -- both on rows and columns
        # finally, convert the result back into a triangular matrix.
        return np.triu(np.dot(np.dot(P, full_map), P.T))

    def _calculate_map(self):
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

            try:
                for ir in v[True]:
                    for ic in v[False]:
                        if ic > ir:
                            _map[ir][ic] += 1
                        else:
                            _map[ic][ir] += 1
            except IndexError as er:
                print er
                print ic, ir, v
                sys.exit(1)

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
            start = int(self.offsets['delta'][i] / self.bin_width)
            if i < len(self.offsets) - 1:
                end = int(self.offsets['delta'][i + 1] / self.bin_width)
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

    @staticmethod
    def plot_map(file_name, cmap, bin_count, remove_diag=False):
        cmap = cmap.astype(np.float64)
        # being pure counts, we just add 1 to every element. After log transformation
        # the empty elements will be zero.
        cmap += 1.0

        if remove_diag:
            np.fill_diagonal(cmap, np.min(cmap))

        fig = plt.figure(frameon=False)
        img_h = float(bin_count) / 100
        fig.set_size_inches(img_h, img_h)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(np.log(cmap), interpolation='nearest')
        fig.savefig(file_name, dpi=100)

    @staticmethod
    def write_map(file_name, cmap):
        np.savetxt(file_name, cmap)


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

    parser = argparse.ArgumentParser(description='Create a HiC fragment map from mapped reads')
    parser.add_argument('-N', '--max-iter', type=int, default=1, help='Number of brute forace iterations')
    parser.add_argument('--ref-seq', help='Reference sequence')
    parser.add_argument('--enzyme', help='Enzyme used in HiC restriction digest')
    parser.add_argument('--simu-reads', default=False, action='store_true', help='Handle simulator reads')
    parser.add_argument('--bin-width', type=int, default=10, help='Bin size in bp (25000)')
    parser.add_argument('--remove-diag', default=False, action='store_true',
                        help='Remove the central diagonal from plot')
    parser.add_argument('refseq', metavar='FASTA', help='Reference sequence')
    parser.add_argument('bamfile', metavar='BAMFILE', help='BAM file to read')
    parser.add_argument('output', metavar='OUTPUT_BASE', nargs=1, help='Output base name')
    args = parser.parse_args()

    fm = FragmentMap(args.bamfile, args.refseq, args.enzyme,
                     bin_width=args.bin_width, simu_reads=args.simu_reads)

    print 'Writing raw output'
    fm.write_map('{0}.raw.cm'.format(args.output[0]), fm.raw_map)
    fm.plot_map('{0}.raw.png'.format(args.output[0]),
                fm.raw_map, fm.groupings.total_bins(), remove_diag=args.remove_diag)

    with open(args.output[0] + '.log', 'w') as log_h:
        starting_order = fm.order.order.tolist()
        max_n = args.max_iter
        max_lL = None
        max_order = None
        max_names = None
        for n in xrange(max_n):
            if n % 250 == 0:
                print 'Finished {0}/{1} calcs'.format(n, max_n)

            #logL = fm.calc_likelihood(1.0e-8, 250000.0, 0.11)
            logL = fm.calc_likelihood(2.3e-6, 430000.0, 0.11)
            if logL > max_lL or not max_lL:
                max_lL = logL
                max_order = fm.order.order.copy()
                print max_lL
                print max_order
                max_names = np.array(fm.order.names)[fm.order.order].tolist()

            order_str = ' '.join([str(oi) for oi in fm.order.order])
            log_h.write('{0} {1}\n'.format(n, max_names))
            log_h.write('{0} {1} {2}\n'.format(n, order_str, logL))
            fm.shuffle_order()

        order_str = ' '.join([str(oi) for oi in max_order])
        print 'Maximum logL was {0} for order {1}'.format(max_lL, order_str)

        with open('{0}.reord.fasta'.format(args.output[0]), 'w') as out_h:
            seqs = [seq for seq in SeqIO.parse(args.refseq, 'fasta')]
            for oi in max_order:
                SeqIO.write(seqs[oi], out_h, 'fasta')
