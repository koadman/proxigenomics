#!/usr/bin/env python
import argparse
import pysam
import numpy as np
import networkx as nx
import community as com
import math
import sys
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch
from scipy.misc import factorial
from scipy.stats import poisson
from scipy.stats import geom

# use matplotlib without x-server
#import matplotlib as mpl
#mpl.use('Agg')
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


def inverse_edge_weights(g):
    for u, v in g.edges():
        g.edge[u][v]['weight'] = 1.0/g[u][v]['weight']


def edgeiter_to_nodelist(edge_iter):
    nlist = []
    for ei in edge_iter:
        for ni in ei:
            if ni not in nlist:
                nlist.append(ni)
    return nlist


def dfs_weighted(g, source=None):
    if source is None:
        # produce edges for all components
        nodes = g
    else:
        # produce edges for components with source
        nodes = [source]
    visited = set()
    for start in nodes:
        if start in visited:
            continue
        visited.add(start)
        stack = [(start, iter(sorted(g[start], key=lambda x: -g[start][x]['weight'])))]
        while stack:
            parent, children = stack[-1]
            try:
                child = next(children)
                if child not in visited:
                    yield parent, child
                    visited.add(child)
                    stack.append((child, iter(sorted(g[child], key=lambda x: -g[child][x]['weight']))))
            except StopIteration:
                stack.pop()

def decompose_graph(g):
    """
    Using the Louvain algorithm for community detection, as
    implemented in the community module, determine the partitioning
    which maximises the modularity. For each individual partition
    create the sub-graph of g

    :param g: the graph to decompose
    :return: the set of sub-graphs which form the best partitioning of g
    """
    decomposed = []
    part = com.best_partition(g)
    part_labels = np.unique(part.values())

    # for each partition, create the sub-graph
    for pi in part_labels:
        # start with a complete copy of the graph
        gi = g.copy()
        # build the list of nodes not in this partition and remove them
        to_remove = [n for n in g.nodes_iter() if part[n] != pi]
        gi.remove_nodes_from(to_remove)
        decomposed.append(gi)

    return decomposed

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
    def digestion_sites(seq_list, enzyme_names=[], min_sites=1):
        """
        Return a list of sites per sequence, preserving the input list order.
        :param seq_list: list of sequences to analyze
        :param enzyme_names: enzyme used in digestion
        :param min_sites: minimum sites required for a sequence to be included.
        :param min_length: minimum sequence length to be included.
        :return: list of sites per sequence
        """
        sites = []
        for seq in seq_list:
            if seq['excluded']:
                continue

            ds = DigestedSequence(enzyme_names, seq['record'])
            seq_sites = ds.get_sites()

            if len(seq_sites) < min_sites:
                print '\tExcluded {0} (length {1}) with only {2} sites'.format(seq['record'].id, len(seq['record']), len(seq_sites))
                #continue

            sites.append({'name': seq['record'].id, 'pos': seq_sites})
        return sites


class Grouping:

    # masked bit
    MASK = -1

    def __init__(self, site_list, bin_width, min_sites):
        total_sites = 0
        self.bins = []
        self.map = []
        self.bin_width = bin_width
        self.min_sites = min_sites

        for seq_sites in site_list:

            max_n = len(seq_sites['pos'])

            if max_n == 0:
                # Sequences which had no sites get empty place holders
                # this maintains order of lists
                self.bins.append(Grouping.MASK)
                self.map.append(None)
                continue

            # required bins is the quotient, but the least number
            # of bins is 1.
            required_bins = max_n / bin_width if max_n >= bin_width else 1
            if max_n > bin_width and float(max_n % bin_width)/bin_width >= 0.5:
                # add an extra bin when there is a sufficient number of left-overs
                # but do not split a single bin sequence.
                required_bins += 1

            if required_bins == 1 and max_n < min_sites:
                self.bins.append(Grouping.MASK)
                self.map.append(None)
                print '\tSequence had only {0} sites, which is less ' \
                      'than the minimum threshold {1}'.format(required_bins, min_sites)
                continue


            # starting with integer indices, apply a scale factor to adjust
            # them matching the desired binning. This can then be digitised
            # to determine bin membership
            sclidx = np.arange(max_n, dtype=np.int32) * float(required_bins)/max_n

            # determine bin membership, but subtract one from each to get 0-based
            grp = np.digitize(sclidx, np.arange(required_bins)) -1

            self.bins.append(required_bins)
            self.map.append(np.column_stack((seq_sites['pos'], grp)))

            total_sites += max_n

            print '\t{0} sites in {1} bins. Bins: {2}'.format(max_n, required_bins, np.bincount(grp))

        self.bins = np.array(self.bins)

        self.total = total_sites

        # calculate centers of each fragment groupings. This will be used
        # when measuring separation.
        self.centers = []
        for n, num_bins in enumerate(self.bins):
            if num_bins == Grouping.MASK:
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

    def __init__(self, seq_list, site_list):
        """
        Order is initially determined by the order of supplied sequence list. Indices of sequences
        which possessed no restriction sites are indicated in "no_sites".

        :param seq_list: list of sequences
        :param site_list: list of re-sites per sequence.
        """
        self.names = [seq['record'].id for seq in seq_list if not seq['excluded']]
        self.lengths = np.array([len(seq['record']) for seq in seq_list if not seq['excluded']])
        self.order = np.arange(len(self.names))
        self.no_sites = np.array([len(site['pos']) == 0 for site in site_list])

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

    @property
    def active_seq(self):
        return len(self.order.order)

    @property
    def active_len(self):
        return np.sum(self.order.lengths)

    def __init__(self, bam_file, fasta_file, enzyme, min_sites=1, bin_width=3, simu_reads=False):

        min_length = 1000
        self.bin_width = bin_width
        self.min_sites = min_sites
        self.simu_reads = simu_reads

        with open(fasta_file, 'r') as fasta_h:
            seq_list = []
            for seq in SeqIO.parse(fasta_h, 'fasta'):
                excluded = False
                if len(seq) < min_length:
                    print '\tExcluding {0} with short sequence {1} bp'.format(seq.id, len(seq))
                    excluded = True
                seq_list.append({'record': seq, 'excluded': excluded})

            print 'Digesting supplied sequences with {0} ...'.format(enzyme)
            self.sites = DigestedSequence.digestion_sites(seq_list, [enzyme], 5)

            print 'Initializing sequence order ...'
            self.order = SeqOrder(seq_list, self.sites)

        with pysam.AlignmentFile(bam_file, 'rb') as self.bam:
            print 'Analyzing BAM file ...'
            #self.active_seq = len(self.bam.references)
            #self.active_len = sum(self.bam.lengths)
            self.total_reads = self.bam.count()

            print '\tMap based upon mapping containing:\n' \
                  '\t{0} sequences\n' \
                  '\t{1}bp total length\n' \
                  '\t{2} mapped reads'.format(self.active_seq, self.active_len, self.total_reads)

            print 'Binning restriction sites into groups of approximately {0} ...'.format(bin_width)
            self.groupings = Grouping(self.sites, self.bin_width, self.min_sites)

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

    def create_contig_graph(self):
        """
        Create a graph where contigs are nodes and edges linking
        nodes are weighted by the cumulative weight of contacts shared
        between them, normalized by the product of the number of fragments
        involved.

        :return: graph of contigs
        """
        _order = self.order.order
        _bins = self.groupings.bins

        g = nx.Graph()
        g.add_nodes_from(_order)

        for i in xrange(len(_order)):
            if _bins[i] == Grouping.MASK:
                continue
            for j in xrange(i+1, len(_order)):
                if _bins[j] == Grouping.MASK:
                    continue
                sm = self.get_submatrix(i, j) + 1
                n_frg = len(self.groupings.map[i]) * len(self.groupings.map[j])
                # networkx written to graphml will choke on numpy floats
                w = float(np.sum(sm)) / n_frg
                g.add_edge(i, j, weight=w)
        return g

    def order_contigs_by_hc(self):
        from scipy.cluster.hierarchy import complete
        from scipy.spatial.distance import squareform
        from scipy.cluster.hierarchy import dendrogram
        g = self.create_contig_graph()
        inverse_edge_weights(g)
        D = squareform(nx.adjacency_matrix(g).todense())
        Z = complete(D)
        return dendrogram(Z)['leaves']

    def order_contigs(self):
        """
        Attempt to determine an initial starting order of contigs based
        only upon the cross terms (linking contacts) between each using
        graphical techniques.

        Beginning with a graph of contigs, where edges are weighted by
        contact weight, it is decomposed using Louvain modularity. Taking
        inverse edge weights, the shortest path of the minimum spanning
        tree of each subgraph is used to define an order. The subgraph
        orderings are then concatenated together to define a full
        ordering of the sample.

        Those with no edges, are included by appear in an indeterminate
        order.

        :return: order of contigs
        """
        g = self.create_contig_graph()
        decomposed_subgraphs = decompose_graph(g)

        isolates = []
        new_order = []
        for gi in decomposed_subgraphs:
            if gi.order() > 1:
                inverse_edge_weights(gi)
                mst = nx.minimum_spanning_tree(gi)
                inverse_edge_weights(gi)
                new_order.extend(edgeiter_to_nodelist(dfs_weighted(mst)))
            else:
                isolates.extend(gi.nodes())

        return new_order + isolates

    def get_submatrix(self, i, j):
        """
        Return the submatrix of the raw contact counts pertaining to the interaction
        between sequence i and j.
        :param i: first sequence id
        :param j: second sequence id
        :return: submatrix (nparray) containing raw counts between these sequences
        """
        _bins = self.groupings.bins
        oi = np.sum(_bins[0:i])
        oj = np.sum(_bins[0:j])
        return self.raw_map[oi:oi+_bins[i], oj:oj+_bins[j]]

    # def calc_likelihood_weight(self, i, j):
    #     if i == j:
    #         raise RuntimeError('Indices cannot be the same. {0},{1}'.format(i,j))
    #
    #     Nd = np.sum(self.raw_map)
    #     L =

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

        # prepare a mask of sequences to skip in calculation
        seq_masked = np.empty(len(self.order.names), dtype=np.bool)
        seq_masked.fill(False)
        seq_masked[self.order.no_sites | self.groupings.bins == Grouping.MASK] = True
        print '{0} sequences will be masked in likelihood calculation'.format(np.sum(seq_masked))

        Nd = np.sum(self.raw_map)
        sumL = 0.0

        num_included_seq = len(self.order.names)
        for i in xrange(num_included_seq):

            #if self.order.no_sites[i]:
            if seq_masked[i]:
                continue

            for j in xrange(i+1, num_included_seq):

                #if self.order.no_sites[j]:
                if seq_masked[j]:
                    continue

                # inter-contig separation defined by cumulative
                # intervening contig length.
                L = self.order.intervening(i, j)

                # bin centers
                centers_i = self.groupings.centers[i]
                centers_j = self.groupings.centers[j]

                # determine relative origin for measuring separation
                # between sequences. If i comes before j, then distances
                # to j will be measured from the end of i -- and visa versa
                if self.order.is_first(i, j):
                    s_i = self.order.lengths[i] - centers_i
                    s_j = centers_j
                else:
                    s_i = centers_i
                    s_j = self.order.lengths[j] - centers_j

                #
                d_ij = np.abs(L + s_i[:, np.newaxis] - s_j)

                # Here we are using a peice-wise continuous function defined in GRAAL
                # to relate separation distance to Poisson rate parameter mu.
                # Above a certain separation distance, the probability reaches a constant minimum value.
                q_ij = np.piecewise(d_ij, [d_ij < 3e6, d_ij > 3e6],
                                    [lambda x: 0.5*(1.0 / 3e6 - (1. - 6e-6)**x * np.log(1.0 - 6e-6)),
                                     1.0e-8])

                n_ij = self.get_submatrix(i, j)
                p_ij = poisson.logpmf(n_ij, mu=(Nd * q_ij))
                sumL += np.sum(p_ij)
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

            print '\tFetching reads for {0} ...'.format(seq_name)
            for r in self.bam.fetch(seq_name):
                n += 1
                if n % print_rate == 0:
                    msg = 'Processing {0}/{1} {2}'.format(ctg_idx+1, self.active_seq, seq_name)
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
              'representing {1} bp over {2} sequences'.format(n_bins, self.active_len, self.active_seq)
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
            #print '{0} -> {1} {2} {3},{4} s:{5} {6},{7}'.format(org_i, shuff_i, ointervening_bins, start, stop, shft, curr_bin, curr_bin+_bins[org_i])
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
        if total_pairs < 1000:
            rate = total_pairs
        else:
            rate = total_pairs / 1000

        n = 0
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
            rel_length = self.offsets['length'][i] / float(self.active_len)

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
            #    i+1, self.active_seq, self.offsets['name'][i], start, end, rel_reads, rel_length, scl)
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
    parser.add_argument('--min-sites', type=int, default=1, help='Ignore bins with less than minimum sites [1]')
    parser.add_argument('--enzyme', help='Enzyme used in HiC restriction digest')
    parser.add_argument('--simu-reads', default=False, action='store_true', help='Handle simulator reads')
    parser.add_argument('--bin-width', type=int, default=10, help='Bin size in bp (25000)')
    parser.add_argument('--remove-diag', default=False, action='store_true',
                        help='Remove the central diagonal from plot')
    parser.add_argument('refseq', metavar='FASTA', help='Reference sequence')
    parser.add_argument('bamfile', metavar='BAMFILE', help='BAM file to read')
    parser.add_argument('output', metavar='OUTPUT_BASE', nargs=1, help='Output base name')
    args = parser.parse_args()

    fm = FragmentMap(args.bamfile, args.refseq, args.enzyme, min_sites=args.min_sites,
                     bin_width=args.bin_width, simu_reads=args.simu_reads)

    print 'Writing raw output'
    fm.write_map('{0}.raw.cm'.format(args.output[0]), fm.raw_map)
    fm.plot_map('{0}.raw.png'.format(args.output[0]),
                fm.raw_map, fm.groupings.total_bins(), remove_diag=args.remove_diag)

    with open(args.output[0] + '.log', 'w') as log_h:

        print 'Starting order    ', fm.order.order.tolist()

        init_order = fm.order_contigs()
        print 'Adhoc reordering     ', init_order
        fm.set_order(init_order)
        new_map = fm.reorder_map()
        fm.plot_map('{0}.reorg.png'.format(args.output[0]),
                    new_map, fm.groupings.total_bins(), remove_diag=args.remove_diag)
        logL = fm.calc_likelihood(2.3e-6, 430000.0, 0.11)
        print 'Ad-hoc ordering logL {0}'.format(fm.calc_likelihood(2.3e-6, 430000.0, 0.11))

        hc_order = fm.order_contigs_by_hc()
        print 'HC complete reordering', hc_order
        fm.set_order(hc_order)
        new_map = fm.reorder_map()
        fm.plot_map('{0}.hc.png'.format(args.output[0]),
                    new_map, fm.groupings.total_bins(), remove_diag=args.remove_diag)
        print 'HC ordering logL     {0}'.format(fm.calc_likelihood(2.3e-6, 430000.0, 0.11))

        sys.exit(0)
        max_n = args.max_iter
        max_lL = None
        max_order = None
        max_names = None
        for n in xrange(max_n):
            if n % 50 == 0:
                print 'Finished {0}/{1} calcs'.format(n, max_n)

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

        fm.set_order(max_order)
        new_map = fm.reorder_map()
        fm.plot_map('{0}.maxLL.png'.format(args.output[0]),
                    new_map, fm.groupings.total_bins(), remove_diag=args.remove_diag)
