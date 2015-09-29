#!/usr/bin/env python

from collections import OrderedDict
import sys
import operator
import argparse
import contextlib

import pandas
import networkx


def write_metis(G, metis_file, nodemap_file):
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
    # original node names for downstream tools
    with open(nodemap_file, 'w') as h_map, open(metis_file, 'w') as h_metis:

        # Write header
        h_metis.write('{0} {1} {2}\n'.format(G.number_of_nodes(), G.number_of_edges(), format_type))

        # Write body - the edge connectivity for each node
        for nname, nid in node_map.iteritems():

            # for all adjacent nodes, build the string of connected nodes with weights
            line = []
            for adj, adj_dat in G[nname].iteritems():
                line.extend([str(node_map[adj]), str(int(adj_dat['weight']))])

            h_metis.write(' '.join(line))
            h_metis.write('\n')

            h_map.write('{idx} {name}\n'.format(idx=nid, name=nname))


def check_weights(edgeTable):
    """
    Check that weights in edge table are representable by integers without loss of effective
    precision. Throws a RuntimeError if the conditions are unacceptable.
    :param edgeTable:
    """

    w_max = edgeTable.RAWWEIGHT.max()
    w_min = edgeTable.RAWWEIGHT.min()

    if w_max < 1.0:
        raise RuntimeError('Warning: Metis file format only supports integer weightings. Maximum weight is less than 1')

    elif w_min > 0.0 and w_min < 1.0:
        raise RuntimeError('Warning: Metis file format only supports integer weightings. Minimum weight is between 0 and 1.')


parser = argparse.ArgumentParser(description='Create graph format files from node and edge data')
parser.add_argument('-m', '--minlen', type=int, default=500, help='Minimum contig length')
parser.add_argument('-f', '--fmt', dest='format', default='metis', choices=['metis', 'graphml'],
                    help='Output graph format')
parser.add_argument('--scale-weights', type=float, metavar='FLOAT',
                    help='Scale factor for weights for metis (integer weights only)')
parser.add_argument('edges', metavar='EDGE_CSV', nargs=1, help='Edge csv file')
parser.add_argument('nodes', metavar='NODE_CSV', nargs=1, help='Node csv file')
parser.add_argument('output', metavar='GRAPH_OUT', nargs=1, help='Output file')
parser.add_argument('node_map', metavar='NODEMAP_FILE', nargs='?', help='Output node-map file')
args = parser.parse_args()

try:
    if args.format == 'metis' and args.node_map is None:
        raise RuntimeError('Metis output format requires node_map filename to be specified')

    # Minimum sequence length
    minLength = args.minlen

    # load table of edges
    edgeTable = pandas.read_csv(args.edges[0], sep=' ', dtype={'SOURCE': str, 'TARGET': str, 'RAWWEIGHT': float})
    nodeTable = pandas.read_csv(args.nodes[0], sep=' ', dtype={'ID': str, 'LENGTH': int})

    if args.scale_weights:
        edgeTable.RAWWEIGHT *= args.scale_weights

    for w in edgeTable.RAWWEIGHT:
        if w < 1.0:
            print 'Warning: metis format edge weights must be integer values. {0} will be set to 1'.format(w)
    
    check_weights(edgeTable)

    # upgrade all edges with weights less than 1 to 1.
    edgeTable.loc[edgeTable.RAWWEIGHT < 1.0, 'RAWWEIGHT'] = 1.0
    
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
        write_metis(G, args.output[0], args.node_map)
    elif args.format == 'graphml':
        networkx.write_graphml(G, args.output[0])
    else:
        raise RuntimeError('unsupported graph format requested')

except RuntimeError as er:
    sys.stderr.write('Error: {0}\n'.format(er.message))
    sys.exit(1)
