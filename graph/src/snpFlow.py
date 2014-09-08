import networkx as nx
import numpy as np
import pandas as pd
import sys


class SNPNode:

    def __init__(self, str_rep):
        tok = str_rep.split(':')
        self.contig = tok[0]
        self.position = int(tok[1])

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return '{0}:{1}'.format(self.contig, self.position)

    def __hash__(self):
        return hash((self.contig, self.position))

    # order by contig then position
    def __cmp__(self, other):
        if not isinstance(other, SNPNode):
            return NotImplemented
        if self.contig == other.contig:
            return cmp(self.position, other.position)
        else:
            return cmp(self.contig, other.contig)


class SplitNode(object):

    # 1-based coords, end inclusive
    def __init__(self, contig, start, end=None):
        self.contig = contig
        self.start = start
        if end is None:
            self.end = start
        else:
            self.end = end

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return '{0}:{1}:{2}'.format(self.contig, self.start, self.end)

    def __hash__(self):
        return hash((self.contig, self.start, self.end))

    def __eq__(self, other):
        if type(other) is not type(self):
            return NotImplemented
        return other.contig == self.contig and other.start == self.start and other.end == self.end

    def __ne__(self, other):
        return not self.__eq__(other)


class VariantSplitNode(SplitNode):

    def __init__(self, contig, pos, base):
        super(VariantSplitNode, self).__init__(contig, pos)
        self.base = base

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return '{0}:{1}'.format(super(VariantSplitNode, self).__repr__(), self.base)

    def __hash__(self):
        return hash((super(VariantSplitNode, self).__hash__(), self.base))

    def __eq__(self, other):
        if type(other) is not type(self):
            return NotImplemented
        return super(VariantSplitNode, self).__eq__(other) and other.base == self.base

    def __ne__(self, other):
        return not self.__eq__(other)


if len(sys.argv) != 3:
    print 'Usage: [snp graph] [graph out]'
    sys.exit(1)

snp_graph = nx.read_graphml(sys.argv[1], node_type=SNPNode)

split_graph = nx.Graph()
for e in snp_graph.edges_iter():

    # original source
    u = SNPNode(str(e[0]))

    # new split source (ref and var)
    udat = snp_graph.node[e[0]]
    uR = VariantSplitNode(u.contig, u.position, udat['ref'])
    uV = VariantSplitNode(u.contig, u.position, udat['var'])

    # original target
    v = SNPNode(str(e[1]))

    # new split target (ref and var)
    vdat = snp_graph.node[e[1]]
    vR = VariantSplitNode(v.contig, v.position, vdat['ref'])
    vV = VariantSplitNode(v.contig, v.position, vdat['var'])

    # weights for each 4 way are a combination of ratios between Ref and Var.
    split_graph.add_edge(uR, vR, weight=(1-udat['ratio']) * (1-vdat['ratio']))
    split_graph.add_edge(uV, vR, weight=udat['ratio'] * (1-vdat['ratio']))
    split_graph.add_edge(uR, vV, weight=(1-udat['ratio']) * vdat['ratio'])
    split_graph.add_edge(uV, vV, weight=udat['ratio'] * vdat['ratio'])

# For each connected component of the graph
#   calculate the pairwise maximum flow where capacity is set equal to snp edge weights.
#   write each pairwise set as a matrix.
#
for n, sg in enumerate(nx.connected_component_subgraphs(split_graph, copy=True)):
    flows = None
    header = []
    for i, u in enumerate(sg.nodes_iter()):
        header.append(str(u))
        fr = []
        for j, v in enumerate(sg.nodes_iter()):
            if i >= j:
                # calculating only half of the symmetric matrix
                fr.append(0)
                continue
            fv, fd = nx.maximum_flow(sg, u, v, capacity='weight')
            fr.append(fv)

        if flows is None:
            flows = np.array(fr, dtype=float)
        else:
            flows = np.vstack((flows, fr))

    # to get a sensibly formatted matrix, use Pandas
    df = pd.DataFrame(flows, columns=header, index=header)
    df.to_csv('sg{0}.txt'.format(n), sep=' ', float_format='%.4f')

nx.write_graphml(split_graph, sys.argv[2])
