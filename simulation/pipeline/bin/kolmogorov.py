#!/usr/bin/env python

import networkx as nx
import argparse
import zlib

def kolmogorov(s):
    """
    Approximate a measure of Kolmogorov complexity by using lossless compression
    to compare input and output of a string. Here, a graph can be represented
    as a linearised adjacency matrix.

    :param s: a string
    :return: measure of complexity as defined by the degree to which a string can be compressed without loss of
    information [0..1].
    """
    if len(s) < 100:
        raise RuntimeError('Strings shorter than 100 bytes are not effectively '
                           'compressed by zlib, therefore we cannot proceed')
    s_compressed = zlib.compress(s, 9)
    return float(len(s_compressed))/float(len(s))

parser = argparse.ArgumentParser(description="Kolmogorov complexity estimate for a graph")
parser.add_argument('input', help='GraphML format graph file to analyse')
args = parser.parse_args()

print '  Reading graph {0}'.format(args.input)
g = nx.read_graphml(args.input)

print '  Determining adjacency matrix ...'
g_adj = nx.adjacency_matrix(g)
g_adj = g_adj.todense()

# convert the adjacency matrix (which might be weighted) into a binary string.
# where 0 no connection, 1 connection with weight > 0.
print '  Converting adjacency matrix to binary string representation ...'
g_str = ''.join(map(str, g_adj.flatten().astype(bool).astype(int).tolist()[0]))


print '  Compressing string ...'
print 'Kolmogorov complexity estimate: {0}'.format(kolmogorov(g_str))
