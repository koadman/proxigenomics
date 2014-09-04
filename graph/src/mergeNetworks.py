import networkx as nx
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


if len(sys.argv) != 4:
    print 'Usage: [scf graph] [snp graph] [graph out]'
    sys.exit(1)

ctg_graph = nx.read_graphml(sys.argv[1])
snp_graph = nx.read_graphml(sys.argv[2], node_type=SNPNode)

split_graph = nx.Graph()
for ctgNode, ctgData in ctg_graph.nodes_iter(data=True):

    # build the list of snps in a single contig (with assocatiated data)
    snp_list = [(n, dat) for n, dat in snp_graph.nodes_iter(data=True) if n.contig == ctgNode]

    if len(snp_list) <= 0:
        continue

    # order by position asc
    snp_list.sort(key=lambda tup: tup[0])

    prevFrg = None
    x0 = 1
    for s in snp_list:

        # each SNP produces 3 fragments
        frgA = SplitNode(ctgNode, x0, s[0].position-1)
        frgRef = VariantSplitNode(ctgNode, s[0].position, s[1]['ref'])
        frgVar = VariantSplitNode(ctgNode, s[0].position, s[1]['var'])

        # 1 segment becomes 3
        split_graph.add_node(frgA)
        split_graph.add_node(frgRef, vtype='ref')
        split_graph.add_node(frgVar, vtype='var')

        # join the nodes - ctg to ref/var
        split_graph.add_edge(frgA, frgRef)
        split_graph.add_edge(frgA, frgVar)

        # attach new 3-node tree to daisy-chain
        if prevFrg is not None:
            split_graph.add_edge(prevFrg[0], frgA)
            split_graph.add_edge(prevFrg[1], frgA)

        prevFrg = [frgRef, frgVar]
        x0 = s[0].position + 1

    # last segment corner-case
    if x0 < ctgData['length']:
        frg = SplitNode(ctgNode, x0, ctgData['length'])
        split_graph.add_node(frg)
        split_graph.add_edge(prevFrg[0], frg)
        split_graph.add_edge(prevFrg[1], frg)

# Reintroduce the SNP-SNP links, where now a single link becomes four.
for e in snp_graph.edges_iter():

    # original source
    u = SNPNode(str(e[0]))

    # new split target (ref and var)
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

nx.write_graphml(split_graph, sys.argv[3])
