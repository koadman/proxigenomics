#!/usr/bin/env python

from collections import OrderedDict
import sys
import operator
import argparse

import pandas
import networkx


def write_metis(G, metis_file):
    """Metis format
    - this format imposes implicit limitation that nodes names must be integers
    - provision for weights also appears to be limited to integers
    """
    # sort nodes by original names and create a set of synthetic numeric identifiers
    # to comply with Metis format limitations
    sorted_nodes = sorted(networkx.degree(G).iteritems(), key=operator.itemgetter(0))
    node_map = OrderedDict()
    for nid, nname in enumerate(sorted_nodes, start=1):
        node_map[nname[0]] = nid

    # Metis format uses a binary flag field in header to specify weights.
    # '1' means just edge weights
    format_type = '1'

    # write out two files.
    # - a list of node names to node ids
    # - a metis format file
    try:
        # original node names for downstream tools
        map_file = open('{0}.nodemap'.format(metis_file), 'w')
        metis_file = open(metis_file, 'w')

        try:
            # Header
            metis_file.write('{0} {1} {2}\n'.format(G.number_of_nodes(), G.number_of_edges(), format_type))

            # Body - the edge connectivity for each node
            for nname, nid in node_map.iteritems():

                # for all adjacent nodes, build the string of connected nodes with weights
                line = []
                for adj, adj_dat in G[nname].iteritems():
                    line.extend([str(node_map[adj]), str(adj_dat['weight'])])

                metis_file.write(' '.join(line))
                metis_file.write('\n')

                map_file.write('{idx} {name}\n'.format(idx=nid, name=nname))
        finally:
            map_file.close()
            metis_file.close()
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

# Minimum sequence length
minLength = args.minlen

# load table of edges
edgeTable = pandas.read_csv(args.edges[0], sep=' ')
nodeTable = pandas.read_csv(args.nodes[0], sep=' ')

# Remove sequences below length threshold
filteredIDs = nodeTable[nodeTable.LENGTH > minLength].ID
filteredEdges = edgeTable[edgeTable.TARGET.isin(filteredIDs) & edgeTable.SOURCE.isin(filteredIDs)]

# for each node, get the list of neighbours and their raw weights
ndict = {}
for n in filteredIDs:
    links = {}
    # the neighbours of node 'n' when listed as TARGETs
    for e in filteredEdges[filteredEdges.SOURCE == n][['TARGET', 'RAWWEIGHT']].itertuples():
        links[e[1]] = {'weight': float(e[2])}
    # the neighbours of node 'n' when listed as SOURCEs
    for e in filteredEdges[filteredEdges.TARGET == n][['SOURCE', 'RAWWEIGHT']].itertuples():
        links[e[1]] = {'weight': float(e[2])}
    ndict[n] = links

# convert this dict of dicts into a graph object
G = networkx.from_dict_of_dicts(ndict)

if args.format == 'metis':
    write_metis(G, args.output[0])
elif args.format == 'graphml':
    networkx.write_graphml(G, args.output[0])
else:
    print 'Error: unsupported graph format requested'
    sys.exit(1)
