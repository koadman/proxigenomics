#!/usr/bin/env python
#
# Convert a SAM file to edge CSV file suitable for importing into Gephi
#
# There is an assumption that read names contain 'fwd' or 'rev' as a suffix.
#
# Eg. frg01fwd or frg9999rev
#
import sys
import math
import pysam

if len(sys.argv) != 5:
    print 'Usage: [HIC2CTG BAM] [WGS2CTG-IDXSTATS] [edges output] [nodes output]'
    sys.exit(1)


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

# Read the sam file and build a linkage map
linkage_map = {}
with pysam.AlignmentFile(sys.argv[1], 'rb') as bam_file:
    iter_bam = bam_file.fetch()
    for mr in iter_bam:
        if mr.reference_id == -1:
            continue
        read = mr.query_name[:-3]
        rdir = mr.query_name[-3:]

        # We depend on read naming from HiC simulator
        # this could be changed if the simulator tool created
        # read names with conventional Illumina FastQ headers.
        if not (rdir == 'fwd' or rdir == 'rev'):
            raise RuntimeError('Reads in alignment file do not conform to expected convention '
                                '[a-zA-Z]+[0-9]+(fwd|ref)')

        contig = bam_file.getrname(mr.reference_id)
        linkage = linkage_map.get(read)
        if linkage is None:
            linkage_map[read] = [(contig, rdir)]
        else:
            linkage.append((contig, rdir))


def update_node_map(line):
    fields = line.rsplit()
    if len(fields) != 4:
        raise Exception('Invalid line in idxstats file [' + line + ']')
    node_map[fields[0]] = Node(fields[0], int(fields[1]), int(fields[2]))

# Read the idxstats file and build the node list
node_map = {}
with open(sys.argv[2], 'r') as h_in:
    [update_node_map(line) for line in h_in if not line.startswith('*')]

# From the set of all linkages, convert this information
# into inter-contig edges, where the nodes are contigs.
# Count the number of redundant links as a raw edge weight.
edge_map = {}
for (insert, linkage) in linkage_map.iteritems():
    for i in range(len(linkage)):
        for j in range(i):
            if linkage[i][1] != linkage[j][1] and linkage[i][0] != linkage[j][0]:
                e = Edge([node_map[linkage[i][0]], node_map[linkage[j][0]]])
                if e.get_id() not in edge_map:
                    edge_map[e.get_id()] = e
                else:
                    e = edge_map.get(e.get_id())
                    e.inc_weight()

# Write out all edges to stdout
# Extra columns for the Gephi
with open(sys.argv[3], 'w') as h_out:
    h_out.write("SOURCE TARGET RAWWEIGHT WEIGHT TYPE\n")
    for e in edge_map.values():
        h_out.write('{0} UNDIRECTED\n'.format(e))

with open(sys.argv[4], 'w') as h_out:
    h_out.write('ID LENGTH READS\n')
    for n in node_map.values():
        h_out.write('{0}\n'.format(n))
