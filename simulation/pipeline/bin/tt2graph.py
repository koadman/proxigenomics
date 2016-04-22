#!/usr/bin/env python
import truthtable as tt
import networkx as nx
import sys
#import itertools

if len(sys.argv) != 3:
    print 'usage: <mcl> <out>'
    sys.exit(1)

t = tt.read_truth(sys.argv[1])
t = t.invert()

g = nx.DiGraph()

for ci, seq_list in t.iteritems():
    g.add_node(ci, type='hub')
    for si in seq_list:
        g.add_node(si, type='child')
        g.add_edge(ci, si)

    

#for clid, seqs in t.iteritems():
#    print clid, len(seqs)
#    for si in seqs:
#        for sj in seqs:
#            g.add_edge(si, sj
    #it = itertools.combinations(seqs, 2)
    #for ei in it:
    #    g.add_edge(*ei, clid=clid)

nx.write_graphml(g,sys.argv[2])

