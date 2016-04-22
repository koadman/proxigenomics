#!/usr/bin/env python
import networkx as nx
import argparse
import math

parser = argparse.ArgumentParser(description='Normalise graph')
parser.add_argument('-b', '--threshold', default=0.0, type=float, help='Edge weight threshold to apply before clustering')
parser.add_argument('input', help='Input graph (graphml format)')
parser.add_argument('output', help='Output file')
args = parser.parse_args()

g = nx.read_graphml(args.input)

thres_g = g.copy()
# thresholding, current normalises first because we've got raw weights
min_w = None
for u,v,d in g.edges_iter(data=True):

    w = float(d['weight'])
    l_u = float(g.node[u]['length']) #/256.
    l_v = float(g.node[v]['length']) #/256.
    harm_mean = (l_u*l_v)/(l_u + l_v)
    thres_g[u][v]['weight'] = w / math.log10(harm_mean)
    
#    thres_g[u][v]['weight'] = 65536.0 * float(d['weight']) / float(g.node[u]['length'] * g.node[v]['length'])

    if thres_g[u][v]['weight'] < min_w or not min_w:
        min_w = thres_g[u][v]['weight']


# bring the lowest weight up to 1. Helps with some tools that numbers
# are not too small
scl = 1.0 / min_w
for u,v,d in thres_g.edges_iter(data=True):
    thres_g[u][v]['weight'] *= scl
    if thres_g[u][v]['weight'] < args.threshold:
        thres_g.remove_edge(u, v)


if args.threshold > 0.0:
    print 'Thresholding reduced graph edge count from {0} to {1}'.format(g.size(), thres_g.size())

nx.write_graphml(thres_g, args.output)

