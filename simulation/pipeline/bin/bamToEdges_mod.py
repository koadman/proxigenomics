#!/usr/bin/env python
#
# Convert a SAM file to edge CSV file suitable for importing into Gephi
#
# There is an assumption that read names contain 'fwd' or 'rev' as a suffix.
#
# Eg. frg01fwd or frg9999rev
#
import math
import pysam
import argparse
import re
import sys
import networkx as nx


class Edge:
    """Represents an edge in the network of contigs linked
    by Hi-C read pairs.
    """

    def __init__(self, nodes=[]):
        nodes.sort()
        self.nodes = nodes
        self.weight = 1

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.nodes[0] == other.nodes[0] and self.nodes[1] == other.nodes[1]

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return '{source} {target} {weight} {normweight}'.format(
                source=self.nodes[0].id,
                target=self.nodes[1].id,
                weight=str(self.weight),
                normweight=str(self.norm_weight()))

    def inc_weight(self):
        self.weight += 1

    def get_id(self):
        return self.nodes[0].id + self.nodes[1].id

    def norm_weight(self):
        return self.weight / math.sqrt(self.nodes[0].length * self.nodes[1].length)


class Node:
    """Represents a node in the network
    """

    def __init__(self, id, length, reads):
        self.id = id
        self.length = length
        self.reads = reads

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.id == other.id

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self.id < other.id

    def __gt__(self, other):
        return self.id > other.id

    def __le__(self, other):
        return self.id <= other.id

    def __ge__(self, other):
        return self.id >= other.id

    def __str__(self):
        return '{id} {length} {reads}'.format(
            id=self.id, length=self.length, reads=self.reads)


def update_linkage_map(l):
    """Parse the line for new information about contig linkages. These
    may be self-self linkages or between inter-contig.
    """
    field = l.rstrip('\n').lstrip().split()
    read = field[0][:-3]
    rdir = field[0][-3:]
    contig = field[2]
    linkage = linkage_map.get(read)
    if linkage is None:
        linkage_map[read] = [(contig, rdir)]
    else:
        linkage.append((contig, rdir))


# Filter lines beginning with '@' and any line where the
# subject sequence is listed as '*'.
def filter(line):
    if line.startswith('@'): return True
    fields = line.rsplit()
    if fields[2] == '*': return True
    return False

def split_name(query_name):
    """
    Following our HiC naming convention, split query names into their
    components: (id, direction)

    :param query_name: query name to split
    :return: (id, direction)
    """
    return query_name[:-3], query_name[-3:]


cig_matcher = re.compile(r'([0-9]+)M')
def count_matches(cigar):
    """
    Sum the length of aligned regions listed in CIGAR pattern
    :param cigar: SAMtools CIGAR string
    :return: total number of aligned bases
    """
    return sum([int(seg_len) for seg_len in cig_matcher.findall(cigar)])


NOT_ALLOWED = set([2,3,6,8])
def strong_match(mr, min_match=None, match_start=True, min_mapq=60):
    """
    Check that a mapped read does not contain disagreeing sequence. Only
    clipped regions are permitted.
    """
    if mr.is_secondary or mr.is_supplementary:
        return False
    
    if len(NOT_ALLOWED & set([t[0] for t in mr.cigartuples])) == 0:
        if match_start and mr.cigartuples[0][0] != 0:
            return False

        if mr.mapping_quality < min_mapq:
            return False

        if min_match:
            n_matches = sum([t[1] for t in mr.cigartuples if t[0] == 0])
            if n_matches < min_match:
                return False

        return True

    else:
        return False


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create edge and node tables from a HiC bam file')
    parser.add_argument('--recover-alts', action='store_true', default=False,
                        help='Recover the alternate alignments from BAM')
    parser.add_argument('--sim', default=False, action='store_true', help='Hi-C simulation read names')
    parser.add_argument('--afmt', choices=['bam', 'psl'], default='bam', help='Alignment file format (bam)')
    parser.add_argument('--minid', type=float, required=False, default=95.0,
                        help='Minimum percentage identity for alignment (95)')
    parser.add_argument('--minlen', type=int, required=False, default=1000,
                        help='Minimum length in bp (1000)')
    parser.add_argument('--mincov', type=float, required=False, default=0.5,
                        help='Minimum coverage of query by alignment (0.5)')

    #parser.add_argument('--wgs', dest='wgs2ctg', metavar='WGS_BAM', nargs=1, help='WGS reads to contigs bam file')
    parser.add_argument('-s', '--add-selfloops', action='store_true', default=False,
                        help='Add self-loops to nodes')
    parser.add_argument('--mapq', default=60, type=int, help='Minimum mapping quality [60]')
    parser.add_argument('--strong', default=None, type=int, help='Accept only mapped reads with no disagreements only clipping')
    parser.add_argument('--graphml', nargs=1, help='Write graphml file')
    parser.add_argument('--split', metavar='BAM', nargs=2, help='Split R1/R2 HiC to contigs bam files')
    parser.add_argument('--merged', metavar='BAM', nargs=1, help='Single merged HiC to contigs bam file')
    parser.add_argument('edge_csv', metavar='EDGE_CSV', nargs=1, help='Edges csv output file')
    parser.add_argument('node_csv', metavar='NODE_CSV', nargs=1, help='Nodes csv output file')
    args = parser.parse_args()

    # TODO We dont need to reference this file.
    # TODO If a scaffold does appear in HiC data then it is not correct to introduce it.

    g = nx.Graph()
    linkage_map = {}

    # input is a BAM file
    if args.afmt == 'bam':

        if args.merged and args.split:
            print 'Error: you must provide either merged OR split bam files'
            sys.exit(1)

        elif args.split and args.split[0] == args.split[1]:
            print 'Error: split files must differ'
            sys.exit(1)

        bam_files = args.split if args.split else args.merged

        for rdir, fn in enumerate(bam_files, start=1):
            print 'Parsing {0}...'.format(fn)

            #rdir = 'fwd' if R == 1 else 'rev'

            with pysam.AlignmentFile(fn, 'rb') as bf:
                for i in xrange(len(bf.references)):
                    rn = bf.references[i]
                    if rn not in g:
                        g.add_node(rn, length=bf.lengths[i])

            reject = 0
            accept = 0
            # Read the sam file and build a linkage map
            with pysam.AlignmentFile(fn, 'rb') as bf:
                iter_bam = bf.fetch()
                for mr in iter_bam:

                    #contig_set = set()

                    if mr.reference_id == -1:
                        reject += 1
                        continue

                    # apply the constraints on alignment coverage and percent identity
                    if args.strong and not strong_match(mr, args.strong, True, args.mapq):
                        reject += 1
                        continue
                    else:
                        perid = dict(mr.cigartuples)[0]/float(mr.query_length) * 100.0
                        cov = mr.query_alignment_length/float(mr.query_length)
                        if cov < args.mincov or perid < args.minid:
                            reject += 1
                            continue

                    accept += 1

                    if args.sim:
                        read = mr.query_name[:-3]
                        rdir = 0  if mr.query_name[-3:] == 'fwd' else 1
                    else:
                        read = mr.query_name

                    # assumed naming convention of HiC simulated reads.
                    #read = mr.query_name[:-3]
                    #rdir = mr.query_name[-3:]
                    #read, rdir = split_name(mr.query_name)
                    #read = mr.query_name

                    # We depend on read naming from HiC simulator
                    # this could be changed if the simulator tool created
                    # read names with conventional Illumina FastQ headers.
                    #if not (rdir == 'fwd' or rdir == 'rev'):
                    #    raise RuntimeError('Reads in alignment file do not conform to expected convention '
                    #                        '[a-zA-Z]+[0-9]+(fwd|ref)')

                    # contig that this alignment line refers to directly
                    #contig_set.add(bf.getrname(mr.reference_id))
                    ctg_name = bf.getrname(mr.reference_id)

                    # if requested, add also alternate alignment contigs for linkage map
                    #if args.recover_alts:
                    #    try:
                            # XA field contains alternate alignments for read, semi-colon delimited
                    #        alts_field = mr.get_tag('XA')
                    #        hit_list = alts_field.split(';')
                            # for each, get the contig name (first element of comma-sep list)
                    #        for hit in hit_list:
                    #            hrec = hit.split(',')
                    #            if len(hrec) != 4:
                    #                continue
                    #            perid = count_matches(hrec[2])/float(mr.query_length) * 100.0
                    #            if perid < args.minid:
                    #                continue
                    #            contig_set.add(hrec[0])
                    #    except KeyError:
                    #        pass

                    #ctg_assocs = [(ctg, rdir) for ctg in contig_set]

                    linkage = linkage_map.get(read)
                    if linkage is None:
                        linkage_map[read] = [(ctg_name, rdir)] #ctg_assocs
                    else:
                        linkage.append((ctg_name, rdir)) #ctg_assocs)
                print 'For {0} -- rejected {1} accepted {2}, rejection rate={3:.1f}%'.format(fn, reject, accept, float(reject)/(reject+accept)*100.)

        print 'Overall -- rejected {0} accepted {1}, rejection rate={2:.1f}%'.format(reject,accept,float(reject)/(reject+accept)*100.) 


    # input is a PSL file
    elif args.afmt == 'psl':

        if args.recover_alts:
            print 'Recovering alternate alignments is only applicable to BAM file parsing'
            sys.exit(1)

        # line marks a non-header line
        psl_dataline = re.compile(r'^[0-9]+\t')

        with open(args.hic2ctg[0], 'r') as h_in:
            for line in h_in:

                # skip header fields
                if not psl_dataline.match(line):
                    continue

                fields = line.rsplit()

                #qname = fields[9]
                # assumed naming convention of HiC simulated reads.
                read, rdir = split_name(fields[9])

                contig = fields[13]
                alen = int(fields[12]) - int(fields[11]) + 1
                qlen = int(fields[10])

                # Taken from BLAT perl script for calculating percentage identity
                matches = int(fields[0])
                mismatches = int(fields[1])
                repmatches = int(fields[2])
                q_num_insert = int(fields[4])
                perid = (1.0 - float(mismatches + q_num_insert) / float(matches + mismatches + repmatches)) * 100.0

                if not g.has_node(contig):
                    g.add_node(contig, length=int(fields[14]))

                # ignore alignment records which fall below mincov or minid
                # wrt the length of the alignment vs query sequence.
                if float(alen)/float(qlen) < args.mincov or perid < args.minid:
                    continue

                linkage = linkage_map.get(read)
                if linkage is None:
                    linkage_map[read] = [(contig, rdir)]
                else:
                    linkage.append((contig, rdir))

    # From the set of all linkages, convert this information
    # into inter-contig edges, where the nodes are contigs.
    # Count the number of redundant links as a raw edge weight.
    edge_map = {}
    unpaired = 0
    odd = 0
    for insert, linkages in linkage_map.iteritems():
        n_links = len(linkages)
        
        if n_links < 2:
            unpaired += 1
            continue
            
        elif n_links > 2:
            odd += 1
            continue
            
        u = linkages[0][0]
        v = linkages[1][0]

        if g.has_edge(u, v):
            g[u][v]['weight'] += 1
        else:
            g.add_edge(u, v, weight=1)
    g.remove_edges_from(g.selfloop_edges())
    print 'unpaired={0} odd={1} order={2} size={3}'.format(unpaired,odd,g.order(),g.size())

#    for insert, linkages in linkage_map.iteritems():
#        for i in range(len(linkages)):
#            for j in range(i):
#
#                # networkx based construction
#                u, ud = linkages[i]
#                v, vd = linkages[j]
#                if g.has_edge(u, v):
#                    g[u][v]['weight'] += 1
#                else:
#                    g.add_edge(u, v, weight=1)

    # insure that there is a self-loop to each node.
    if args.add_selfloops:
        for v in g.nodes():
            if not g.has_edge(v, v):
                g.add_edge(v, v, weight=1)

    # filter nodes on length
    nlist = g.nodes()
    for n in nlist:
        if g.node[n]['length'] < args.minlen:
            g.remove_node(n)

    if args.graphml is not None:
        nx.write_graphml(g, args.graphml[0])

    with open(args.edge_csv[0], 'w') as h_out:
        h_out.write("SOURCE TARGET RAWWEIGHT TYPE\n")
        for u, v, dat in g.edges(data=True):
            h_out.write('{0} {1} {2} UNDIRECTED\n'.format(u, v, dat['weight']))

    with open(args.node_csv[0], 'w') as h_out:
        h_out.write('ID LENGTH\n')
        for v, dat in g.nodes(data=True):
            h_out.write('{0} {1[length]}\n'.format(v, dat))

