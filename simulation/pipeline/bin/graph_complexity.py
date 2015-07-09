#!/usr/bin/env python

import numpy as np
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

def eigen_entropy(g):
    print '  Calculating normalized laplacian ...'
    lap = nx.normalized_laplacian_matrix(g)
    print '  Calculating eigenvalues of laplacian ...'
    eigen_vals = np.linalg.eigvals(lap.A)

    # we will take only non-zero values
    nz_ix = eigen_vals != 0
    if np.sum(nz_ix) <= 0:
        raise RuntimeError('Solving laplacian for eigenvalues returned no non-zero values')
    # entropy calculated on magnitude only
    nz_eig = np.abs(eigen_vals[nz_ix])
    sum_eig = np.sum(nz_eig)
    return -np.sum(nz_eig/sum_eig * np.log(nz_eig / sum_eig))

parser = argparse.ArgumentParser(description='Graph complexity estimation')
parser.add_argument('-m', '--method', choices=['kolmo', 'eigen'], required=True,
                    help='Method to apply')
parser.add_argument('input', help='GraphML format graph file to analyse')
args = parser.parse_args()

print '  Reading graph {0}'.format(args.input)
g = nx.read_graphml(args.input)
print '  Graph contained {0} nodes and {1} edges'.format(g.order(), g.size())

if args.method == 'kolmo':
    print '  Determining adjacency matrix ...'
    g_adj = nx.adjacency_matrix(g)
    g_adj = g_adj.todense()

    # convert the adjacency matrix (which might be weighted) into a binary string.
    # where 0 no connection, 1 connection with weight > 0.
    print '  Converting adjacency matrix to binary string representation ...'
    g_str = ''.join(map(str, g_adj.flatten().astype(bool).astype(int).tolist()[0]))

    try:
        k = kolmogorov(g_str)
    except RuntimeError as er:
        print er
        k = None
    print 'Kolmogorov complexity estimate: {0}'.format(k)

elif args.method == 'eigen':
    try:
        S = eigen_entropy(g)
    except RuntimeError as er:
        print er
        S = None
    print 'Eigen value entropy: {0}'.format(S)
