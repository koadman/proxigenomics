#!/usr/bin/env python

from collections import OrderedDict
import sys
import operator
import argparse

import pandas
import networkx


def write_metis(G, graph_out):
    """Metis format
    - this format imposes implicit limitation that nodes names must be integers
    - provision for weights also appears to be limited to integers
    """
    # sort nodes by original names
    nsrt = sorted(networkx.degree(G).iteritems(), key=operator.itemgetter(0))

    # assign integer ids to all nodes
    node_map = OrderedDict()
    for i, v in enumerate(nsrt):
        node_map[v[0]] = i + 1

    format_type = '1'  # edge weights only

    # write out two files.
    # - a list of node names to node ids
    # - a metis format file
    try:
        map_out = open('{0}.nodemap'.format(graph_out), 'w')
        graph_out = open(graph_out, 'w')

        try:
            # header to metis is "#nodes #edges fmt", where fnt is a format code
            # specifying whether or not the graph is weighted (unweighted, node weights,
            # edge weights or both)
            graph_out.write(
                '{nodes} {edges} {fmt}\n'.format(nodes=G.number_of_nodes(), edges=G.number_of_edges(), fmt=format_type))
            for k, v in node_map.iteritems():
                line = []
                for l, m in G[k].iteritems():
                    line += [str(node_map[l]), str(m['weight'])]
                graph_out.write(' '.join(line) + '\n')
                map_out.write('{idx} {name}\n'.format(idx=v, name=k))
        finally:
            map_out.close()
            graph_out.close()
    except IOError, ex:
        print "Error:", ex

parser = argparse.ArgumentParser(description='Create graph format files from node and edge data')
parser.add_argument('-m', '--minlen', type=int, default=500, help='Minimum contig length')
parser.add_argument('-f', '--fmt', dest='format', default='metis', choices=['metis', 'graphml'],
                    help='Output graph format')
parser.add_argument('edges', metavar='EDGE_CSV', nargs=1, help='Edge csv file')
parser.add_argument('nodes', metavar='NODE_CSV', nargs=1, help='Node csv file')
parser.add_argument('output', metavar='OUTPUT', nargs=1, help='Output file')
args = parser.parse_args()

#if len(sys.argv) != 6:
#    print 'Usage [min length] [edge csv] [node csv] [nodemap out] [metis out]'
#    sys.exit(1)

# Minimum sequence length
minLength = args.minlen

# load table of edges
edgeTable = pandas.read_csv(sys.argv[2], sep=' ')
nodeTable = pandas.read_csv(sys.argv[3], sep=' ')

# Remove sequences below length threshold
filteredIDs = nodeTable[nodeTable.LENGTH > minLength].ID
filteredEdges = edgeTable[edgeTable.TARGET.isin(filteredIDs) & edgeTable.SOURCE.isin(filteredIDs)]

# for each node, get the list of neighbours and their raw weights
ndict = {}
for n in filteredIDs:
    links = {}
    # the neighbours of node 'n' when listed as TARGETs
    for e in filteredEdges[filteredEdges.SOURCE == n][['TARGET', 'RAWWEIGHT']].itertuples():
        links[e[1]] = {'weight': e[2]}
    # the neighbours of node 'n' when listed as SOURCEs
    for e in filteredEdges[filteredEdges.TARGET == n][['SOURCE', 'RAWWEIGHT']].itertuples():
        links[e[1]] = {'weight': e[2]}
    ndict[n] = links

# convert this dict of dicts into a graph object
G = networkx.from_dict_of_dicts(ndict)

if args.format == 'metis':
    write_metis(G, args.output)
elif args.format == 'graphml':
    networkx.write_graphml(G, args.output)
else:
    print 'Error: unsupported graph format requested'
    sys.exit(1)