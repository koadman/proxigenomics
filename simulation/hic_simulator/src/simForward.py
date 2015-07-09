#!/usr/bin/env python

from Bio import Alphabet
from Bio import SeqIO

from collections import OrderedDict
from optparse import OptionParser

from Bio.Restriction import *

import numpy
import re
import time
import sys

#
# Globals
#

# restriction enzymes used in reaction
GLOBAL_CUT4 = Restriction.NlaIII
GLOBAL_CUT6_1 = Restriction.ClaI
GLOBAL_CUT6_2 = Restriction.SalI

# Mixed geom/unif model
MIXED_GEOM_PROB = 6.0e-6

# Random State from which to draw numbers
# this is initialized at start time
RANDOM_STATE = None


class EmpiricalDistribution:
    """
    Defining an empirical distribution, we can then use it to draw random
    numbers.
    """
    def __init__(self, shape, length, bins=1000):
        self.shape = shape
        self.length = length
        self.scale = 1.0/length
        self.bins = bins
        self.xsample = numpy.linspace(0, self.length, self.bins, endpoint=True, dtype=numpy.float64)
        self.ysample = self._cdf(self.xsample)

    def _cdf(self, x):
        """
        The CDF is the crucial element of this method. Currently, it is hard-coded as
        a linear combination of geometric and uniform distributions. This could easily
        be changed or the class made generic by passing CDF as a function.
        :param x:
        :return:
        """
        return 0.5 * (1.0 - (1.0 - self.shape)**x + self.scale * x)

    def rand(self):
        """
        Using the inverse CDF method, draw a random number for the distribution. This
        method looks up the nearest value of random value x in a sampled representation
        and then interpolates between bin edges. Edge case for the first and last bin
        is to merely use that bin.

        :return: random value following distribution
        """
        xv = self.ysample
        yv = self.xsample
        x = RANDOM_STATE.uniform()
        ix = numpy.searchsorted(xv, x)
        if ix >= self.bins:
            ix -= 1
        elif ix == 0:
            ix += 1
        # interp value
        return yv[ix-1] + (yv[ix] - yv[ix-1]) * ((x - xv[ix-1]) / (xv[ix] - xv[ix-1]))


def find_restriction_sites(enzyme, seq):
    """For supplied enzyme, find all restriction sites in a given sequence
    returns list of sites.
    """
    return enzyme.search(seq, linear=False)


def find_priming_sites(oligo, seq):
    """For supplied priming sequence, find positions of all matches in a given sequence
    returns list of sites.
    """
    array = []
    for m in re.finditer(oligo, str(seq)):
        array.append(m.end())
    rc_oligo = Seq(oligo)
    rc_oligo.reverse_complement()
    for m in re.finditer(str(rc_oligo), str(seq)):
        array.append(m.end())
    return array


def make_read(seq, fwd_read, length):
    """From sequence, make a forward or reverse read. If the read is longer than the total sequence
    return the entire sequence."""

    # edge case - return whole sequence of read length > sequence
    if length >= len(seq):
        read = seq
        if not fwd_read:
            read = read.reverse_complement(id=True, description=True)
        return read

    # return forward or reverse read
    if fwd_read:
        read = seq[:length]
    else:
        read = seq[-length:]
        read = read.reverse_complement(id=True, description=True)

    return read


def write_reads(handle, sequences, output_format, dummy_q=False):
    """Append sequence to file in fastq format. Currently no effort is spent
    to stream these records"""
    if dummy_q:
        for s in sequences:
            s.letter_annotations['phred_quality'] = [50] * len(s)

    SeqIO.write(sequences, handle, output_format)


class Part:
    """Represents an unligated fragment from one replicon.
    """

    def __init__(self, seq, pos1, pos2, fwd, replicon):
        self.seq = seq
        self.pos1 = pos1
        self.pos2 = pos2
        self.fwd = fwd
        self.replicon = replicon

    def __repr__(self):
        return repr((self.seq, self.pos1, self.pos2, self.fwd, self.replicon))


class Cell:
    """Represents a cell in the community"""

    def __init__(self, name, abundance):
        self.name = name
        self.abundance = float(abundance)
        self.replicon_registry = {}
        self.index_to_name = {}
        self.cdf = None
        self.pdf = None

    def __repr__(self):
        return repr((self.name, self.abundance))

    def __str__(self):
        return self.name

    def init_prob(self):
        """Initialize the probabilities for replicon selection from within a cell. 
        Afterwards, produce a CDF which will be used for random sampling.
        """
        # Number of replicons in cell
        n_rep = len(self.replicon_registry)

        # Uniform probability initially
        prob = numpy.array([1.0 / n_rep] * n_rep)

        # Scale prob by replicon length
        i = 0
        for repA in self.replicon_registry.values():
            self.index_to_name[i] = repA.name
            prob[i] = prob[i] * repA.length()
            i += 1

        # Normalize
        total_prob = sum(prob)
        self.pdf = numpy.divide(prob, total_prob)

        # Initialize the cumulative distribution
        self.cdf = numpy.hstack((0, numpy.cumsum(self.pdf)))

    def register_replicon(self, replicon):
        self.replicon_registry[replicon.name] = replicon

    def number_replicons(self):
        return len(self.replicon_registry)

    def select_replicon(self, x):
        """From a cell, return the index of a replicon by sampling CDF at given value x"""
        idx = numpy.searchsorted(self.cdf, x) - 1
        if idx < 0:  # edge case when sample value is exactly zero.
            return 0
        return idx

    def pick_inter_rep(self, skip_this):
        if self.number_replicons() == 1:
            raise RuntimeError('cannot pick another in single replicon cell')

        # Keep trying until we pick a different rep.
        rep = skip_this
        while rep is skip_this:
            ri = self.select_replicon(RANDOM_STATE.uniform())
            rep = self.replicon_registry.get(self.index_to_name[ri])

        return rep


class Replicon:
    """Represents a replicon which holds a reference to its containing cell"""

    def __init__(self, name, parent_cell, sequence):
        self.name = name
        self.parent_cell = parent_cell
        self.sequence = sequence
        self.cut_sites = {
            '515F': numpy.array(find_priming_sites("GTGCCAGC[AC]GCCGCGGTAA", sequence.seq)),
            '4cut': numpy.array(find_restriction_sites(GLOBAL_CUT4, sequence.seq)),
            '6cut_1': numpy.array(find_restriction_sites(GLOBAL_CUT6_1, sequence.seq)),
            '6cut_2': numpy.array(find_restriction_sites(GLOBAL_CUT6_2, sequence.seq))
        }

        # initialise an empirical distribution for this sequence.
        self.emp_dist = EmpiricalDistribution(MIXED_GEOM_PROB, self.length())

    def __repr__(self):
        return repr((self.name, self.parent_cell, self.sequence))

    def __str__(self):
        return str(self.parent_cell) + '.' + self.name



    def length(self):
        return len(self.sequence)

    def is_alone(self):
        return self.parent_cell.number_replicons() == 1

    def subseq(self, pos1, pos2, fwd):
        """Create a subsequence from replicon where positions are strand relative.
        Those with fwd=False will be reverse complemented.
        """
        # Likely to have a off-by-one error in this code.
        if fwd:
            ss = self.sequence[pos1:pos2]
            ss.description = str(pos1) + "..." + str(pos2) + ":" + str(fwd)
        else:
            ss = self.sequence[pos2:pos1]
            ss.description = str(pos2) + "..." + str(pos1) + ":" + str(fwd)
            ss = ss.reverse_complement(id=True, description=True)
        return ss

    def random_cut_site(self, cut_type):
        """Return a uniformly random cut-site"""
        cs = self.cut_sites[cut_type]
        idx = RANDOM_STATE.randint(low=0, high=len(cs))
        return cs[idx]

    def nearest_cut_site_above(self, cut_type, pos):
        cs = self.cut_sites[cut_type]
        idx = numpy.searchsorted(cs, pos)
        if idx == len(cs):
            return cs[0]
        else:
            return cs[idx]

    def nearest_cut_site_below(self, cut_type, pos):
        cs = self.cut_sites[cut_type]
        idx = numpy.searchsorted(cs, pos)
        if idx == 0:
            return cs[-1]
        else:
            return cs[idx - 1]

    def nearest_cut_site_by_distance(self, cut_type, pos):
        """Find the nearest restriction cut site for the specified cutter type [4cut, 6cut]
        returns genomic position of nearest cut site"""
        cs = self.cut_sites[cut_type]
        idx = numpy.searchsorted(cs, pos)

        # edge cases when insertion point between be before or after
        # beginning or end of linear sequence
        if idx >= len(cs) - 1:
            d1 = (idx - 1, pos - cs[-1])
            d2 = (0, cs[0])
        elif idx == 0:
            d1 = (idx - 1, self.length() - cs[-1])
            d2 = (0, cs[0])
        else:
            d1 = (idx, pos - cs[idx])
            d2 = (idx + 1, cs[idx + 1] - pos)

        if d2[1] < d1[1]:
            return cs[d2[0]]
        else:
            return cs[d1[0]]

    def constrained_upstream_location(self, origin):
        """Return a location (position and strand) on this replicon where the position is
        constrained to follow a prescribed distribution relative to the position of the
        first location.
        
        return location
        """
        delta = self.emp_dist.rand()
        min_len, max_len = 3, self.length()-3
        while delta < min_len or delta > max_len:
            # redraw if too small or large
            delta = self.emp_dist.rand()

        #delta = draw_delta(3, self.length() - 3)
        # TODO The edge cases might have off-by-one errors, does it matter?'
        loc = origin + delta
        if loc > self.length() - 1:
            loc -= self.length()

        return loc


class Community:
    """Represents the community, maintaining a registry of cells and replicons.
    
    Both cells and replicons are expected to be uniquely named across the commmunity. This constraint comes
    from the requirement that a multi-fasta file contain unique sequence identifiers.
    
    A table is used for a basic community definition, with the following three column format. Lines beginning
    with # and empty lines are ignored.
    
    [replicon name] [cell name] [abundance]
    """

    def __init__(self, interrep_prob, table_filename, seq_filename):
        self.pdf = None
        self.cdf = None
        self.totalRawAbundance = 0
        self.replicon_registry = OrderedDict()
        self.index_to_name = None
        self.cell_registry = OrderedDict()
        self.interrep_prob = interrep_prob

        # Read in the sequences
        sequences = SeqIO.to_dict(SeqIO.parse(open(seq_filename), 'fasta', Alphabet.generic_dna))

        # Read table
        with open(table_filename, 'r') as h_table:
            for line in h_table:
                line = line.rstrip().lstrip()
                if line.startswith('#') or len(line) == 0:
                    continue
                field = line.split()
                if len(field) != 3:
                    print 'sequence table has missing fields at [', line, ']'
                    sys.exit(1)
                replicon_name = field[0]
                cell_name = field[1]
                cell_abundance = field[2]
                parent_cell = self.register_cell(cell_name, cell_abundance)
                self.register_replicon(replicon_name, parent_cell, sequences.get(replicon_name))

        # init community wide probs
        self.__init_prob()
        # init individual cell probs
        for cell in self.cell_registry.values():
            cell.init_prob()

    def register_replicon(self, name, parent_cell, sequence):
        """Add a new cell to the community cell registry"""
        replicon = self.replicon_registry.get(name)
        if replicon is None:
            replicon = Replicon(name, parent_cell, sequence)
            self.replicon_registry[name] = replicon
            parent_cell.register_replicon(replicon)
        return replicon

    def register_cell(self, cell_name, cell_abundance):
        """Add a new replicon to the community replicon registry"""
        parent_cell = self.cell_registry.get(cell_name)
        if parent_cell is None:
            parent_cell = Cell(cell_name, cell_abundance)
            self.cell_registry[cell_name] = parent_cell
            self.__update_total_abundance()
        return parent_cell

    # Total abundance of all cells in registry.
    def __update_total_abundance(self):
        """Recalculate the total relative abundance specified for the community by referring to the registry"""
        ab = 0
        for ca in self.cell_registry.values():
            ab += ca.abundance
        self.totalRawAbundance = ab

    def __init_prob(self):
        """Initialize the probabilities for replicon selection, given the abundances, etc.
        Normalization is always applied. Afterwards, produce a CDF which will be used for
        random sampling.
        """
        self.index_to_name = {}
        prob = numpy.zeros(len(self.replicon_registry))
        i = 0
        for repA in self.replicon_registry.values():
            self.index_to_name[i] = repA.name
            prob[i] = repA.parent_cell.abundance / self.totalRawAbundance * repA.length()
            i += 1
        tp = sum(prob)
        self.pdf = numpy.divide(prob, tp)
        'Initialize the cumulative distribution function for the community replicons'
        self.cdf = numpy.hstack((0, numpy.cumsum(self.pdf)))

    def select_replicon(self, x):
        """From the entire community, return the index of a replicon by sampling CDF at given value x"""
        idx = numpy.searchsorted(self.cdf, x) - 1
        if idx < 0:  # edge case when sample value is exactly zero.
            return 0
        return idx

    def pick_replicon(self, skip_index=None):
        """Random selection of replicon from community. If skipIndex supplied, do not return
        the replicon with this index.
        
        return the index"""
        if skip_index is None:
            return self.select_replicon(RANDOM_STATE.uniform())
        else:
            ri = skip_index
            while ri == skip_index:
                ri = self.select_replicon(RANDOM_STATE.uniform())
            return ri

    def is_intrarep_event(self):
        """Choose if the mate is intra or inter replicon associated. This is a simple
        binary paritioning with a chosen threshold frequency.
        """
        return RANDOM_STATE.uniform() > self.interrep_prob

    def get_replicon_by_index(self, index):
        return self.replicon_registry.get(self.index_to_name[index])

    def get_replicon_by_name(self, name):
        return self.replicon_registry.get(name)

    def unconstrained_read_location(self, replicon_index):
        """Return a location (position and strand) on a replicon where the position is
        uniformly sampled across the replicons sequence.
        
        Returns tuple (pos=int, strand=bool)
        """
        replicon = self.get_replicon_by_index(replicon_index)
        return int(RANDOM_STATE.uniform() * replicon.length()), True

    # def constrained_read_location(self, replicon_index, first_location, forward):
    #     """Return a location (position and strand) on a replicon where the position is
    #     constrained to follow a prescribed distribution relative to the position of the
    #     first location.
    #
    #     Strand is always same as first.
    #
    #     return location (pos=int, strand=bool)
    #     """
    #     replicon = self.get_replicon_by_index(replicon_index)
    #     delta = draw_delta(500, replicon.length() - 500)
    #     if forward:
    #         # TODO The edge cases might have off-by-one errors, does it matter?'
    #         loc = first_location + delta
    #         if loc > replicon.length() - 1:
    #             loc -= replicon.length()
    #     else:
    #         loc = first_location - delta
    #         if loc < 0:
    #             loc = replicon.length() - loc
    #     return loc, forward


def make_unconstrained_part_a():
    rep = comm.get_replicon_by_index(comm.pick_replicon())
    pos6c = rep.random_cut_site('6cut_1')
    # pos6c = rep.random_cut_site('515F')
    pos4c = rep.nearest_cut_site_above('4cut', pos6c)
    seq = rep.subseq(pos6c, pos4c, True)
    return Part(seq, pos6c, pos4c, True, rep)


def make_unconstrained_part_b(part_a):
    rep = part_a.replicon.parent_cell.pick_inter_rep(part_a.replicon)
    pos6c = rep.random_cut_site('6cut_2')
    pos4c = rep.nearest_cut_site_below('4cut', pos6c)
    seq = rep.subseq(pos4c, pos6c, True)
    return Part(seq, pos4c, pos6c, True, rep)


def make_constrained_part_b(part_a):
    loc = part_a.replicon.constrained_upstream_location(part_a.pos1)
    pos6c = part_a.replicon.nearest_cut_site_by_distance('6cut_2', loc)
    pos4c = part_a.replicon.nearest_cut_site_below('4cut', pos6c)
    seq = part_a.replicon.subseq(pos4c, pos6c, True)
    return Part(seq, pos4c, pos6c, True, part_a.replicon)


#
# Commandline interface
#
parser = OptionParser()
parser.add_option('-n', '--num-frag', dest='num_frag',
                  help='Number of Hi-C fragments to generate reads', metavar='INT', type='int')
parser.add_option('-l', '--read-length', dest='read_length',
                  help='Length of reads from Hi-C fragments', metavar='INT', type='int')
parser.add_option('-p', '--interrep-prob', dest='inter_prob',
                  help='Probability that a fragment spans two replicons', metavar='FLOAT', type='float')
parser.add_option('-t', '--community-table', dest='comm_table',
                  help='Community profile table', metavar='FILE')
parser.add_option('-s', '--seq', dest='genome_seq',
                  help='Genome sequences for the community', metavar='FILE')
parser.add_option('-r', '--seed', dest='seed',
                  help="Random seed for initialising number generator", metavar='INT', type='int')
parser.add_option('-o', '--output', dest='output_file',
                  help='Output Hi-C reads file', metavar='FILE')
parser.add_option('-f', '--ofmt', dest='output_format', default='fastq',
                  help='Output format', choices=['fasta', 'fastq'], metavar='output_format [fasta, fastq]')
#parser.add_option('--split-reads', dest='split', default=False, action='store_true',
#                  help='Split output reads into separate R1/R2 files')
(options, args) = parser.parse_args()

if options.num_frag is None:
    parser.error('Number of fragments not specified')
if options.read_length is None:
    parser.error('Read length not specified')
if options.inter_prob is None:
    parser.error('Inter-replicon probability not specified')
if options.comm_table is None:
    parser.error('Community profile table not specified')
if options.genome_seq is None:
    parser.error('Genome sequences file not specified')
if options.seed is None:
    options.seed = int(time.time())
if options.output_file is None:
    parser.error('Output file not specified')

#
# Main routine
#

# set state for random number generation
RANDOM_STATE = numpy.random.RandomState(options.seed)

# Initialize community object
print "Initializing community"
comm = Community(options.inter_prob, options.comm_table, options.genome_seq)

# Open output file for writing reads
with open(options.output_file, 'wb') as h_output:

    print "Creating reads"
    skip_count = 0
    overlap_count = 0
    frag_count = 0

    while frag_count < options.num_frag:
        # Fragment creation

        # Create PartA
        # Steps
        # 1) pick a replicon at random
        # 2) pick a 6cutter site at random on replicon
        # 3) flip a coin for strand
        part_a = make_unconstrained_part_a()

        # Create PartB
        # Steps
        # 1) choose if intra or inter replicon
        # 2) if intER create partB as above
        # 3) if intRA select from geometric
        if part_a.replicon.is_alone() or comm.is_intrarep_event():
            # ligation is between two fragments on same replicon
            part_b = make_constrained_part_b(part_a)
        else:
            # ligation crosses replicons
            part_b = make_unconstrained_part_b(part_a)

        # Join parts
        # PartA + PartB
        fragment = part_a.seq + part_b.seq
        if len(fragment) < 200 or len(fragment) > 1000:
            # only accept fragments within a size range
            skip_count += 1
            continue

        if part_b.pos1 < part_a.pos2 and part_a.pos2 < part_b.pos2:
            overlap_count += 1
            continue

        if part_a.pos1 < part_b.pos2 and part_b.pos2 < part_a.pos2:
            overlap_count += 1
            continue

        read1 = make_read(fragment, True, options.read_length)
        read1.id = 'frg{0}fwd'.format(frag_count)
        read1.description = '{0} {1}'.format(part_a.seq.id, part_a.seq.description)

        read2 = make_read(fragment, False, options.read_length)
        read2.id = 'frg{0}rev'.format(frag_count)
        read2.description = '{0} {1}'.format(part_b.seq.id, part_b.seq.description)

        write_reads(h_output, [read1, read2], options.output_format, dummy_q=True)

        frag_count += 1

print "Ignored " + str(skip_count) + " fragments due to length restrictions"
print "Ignored " + str(overlap_count) + " fragments due to overlap"
