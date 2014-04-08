#!/usr/bin/env python
from collections import OrderedDict
import sys
import pandas
import operator
import networkx

def write_metis(G, nodemap_out, metis_out):
	"""Metis format
	- this format imposes implicit limitation that nodes names must be integers
	- provision for weights also appears to be limited to integers
	"""
	# sort nodes by original names
	nsrt = sorted(networkx.degree(G).iteritems(),key=operator.itemgetter(0))

	# assign integer ids to all nodes
	node_map = OrderedDict()
	for i,v in enumerate(nsrt):
		node_map[v[0]] = i+1

	format_type = '1' # edge weights only

	# write out two files.
	# - a list of node names to node ids
	# - a metis format file
	try:
		map_out = open(nodemap_out, 'w')
		graph_out = open(metis_out, 'w')
		try:
			# header to metis is "#nodes #edges fmt", where fnt is a format code
			# specifying whether or not the graph is weighted (unweighted, node weights,
			# edge weights or both)
			graph_out.write('{nodes} {edges} {fmt}\n'.format(nodes=G.number_of_nodes(),edges=G.number_of_edges(),fmt=format_type))
			for k,v in node_map.iteritems():
				line = []
				for l,m in G[k].iteritems():
					line += [str(node_map[l]), str(m['weight'])]
				graph_out.write(' '.join(line) + '\n')
				map_out.write('{idx} {name}\n'.format(idx=v, name=k))
		finally:
			map_out.close()
			graph_out.close()
	except IOError, e:
		print "Error:", e


if len(sys.argv) != 6:
	print 'Usage [min length] [edge csv] [node csv] [nodemap out] [metis out]'
	sys.exit(1)

# Minimum sequence length
minLength = int(sys.argv[1])

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
	for e in filteredEdges[filteredEdges.SOURCE == n][['TARGET','RAWWEIGHT']].itertuples():
		links[e[1]] = {'weight': e[2]}
	# the neighbours of node 'n' when listed as SOURCEs
	for e in filteredEdges[filteredEdges.TARGET == n][['SOURCE','RAWWEIGHT']].itertuples():
		links[e[1]] = {'weight': e[2]}
	ndict[n] = links

# convert this dict of dicts into a graph object
G = networkx.from_dict_of_dicts(ndict)

write_metis(G,sys.argv[4],sys.argv[5])
