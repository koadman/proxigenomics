#!/usr/bin/env python
from Bio import SeqIOfrom Bio.Restriction import *import networkx as nximport numpy as np
import argparseRESTRICTION_BATCH = NoneENZYME = Nonedef count_sites(seq):    sites = ENZYME.search(seqidx[si].seq, linear=True)    return len(sites)parser = argparse.ArgumentParser(description='Filter and normalise graphml file from raw counts')parser.add_argument('--no-self', default=False, action='store_true', help='Remove self-loops')parser.add_argument('-w', '--weight', default=0, type=int, help='Threshold raw edge weight to exclude.[0]')parser.add_argument('-e', '--enzyme', help='Restriction enzyme used')parser.add_argument('cover', help='BBmap coverage file')parser.add_argument('fasta', help='Fasta file for corresponding node sequences')parser.add_argument('graph', help='GraphML format graph file to analyse')parser.add_argument('output', help='Output GraphML file')args = parser.parse_args()g = nx.read_graphml(args.graph)print 'Raw graph contained {0} nodes and {1} edges'.format(g.order(), g.size())
if args.no_self:    g.remove_edges_from(g.selfloop_edges())    print 'Removed self-loops: {0}/{1}'.format(g.order(), g.size())if args.weight > 0:    for u, v, d in g.edges_iter(data=True):        if d['weight'] <= args.weight:            g.remove_edge(u, v)    for u in g.nodes():        if g.degree(u) == 0:            g.remove_node(u)    print 'Removed edges of {0} weight or less: {1}/{2}'.format(args.weight, g.order(), g.size())if args.enzyme:    RESTRICTION_BATCH = RestrictionBatch([args.enzyme])    ENZYME = RestrictionBatch.get(rb, args.enzyme, add=False)seqidx = SeqIO.index(args.fasta, 'fasta')attrib = {}for si in seqidx:    if g.has_node(si):        attr = {}        attr['length'] = len(seqidx[si])        if args.enzyme:            seq_sites = count_sites(seqidx[si])            seq_sites = seq_sites if seq_sites > 0 else 1            attr['site'] = seq_sites        g.node[si].update(attr)with open(args.cover, 'r') as cov_h:    for l in cov_h:        if l.startswith('#ID'):            continue        tok = l.strip().split('\t')        # handle sequences with descriptions        si = tok[0]        if ' ' in si:            si = si.split()[0]        if g.has_node(si):            g.node[si]['depth'] = float(tok[1])            g.node[si]['gc'] = float(tok[8])

for u in g.nodes():
    if g.node[u]['depth'] == 0.0:
        g.remove_node(u)
print 'Removed nodes with zero depth: {0}/{1}'.format(g.order(), g.size())

max_L = float(np.max(np.array(nx.get_node_attributes(g, 'length').values())))
med_D = float(np.median(np.array(nx.get_node_attributes(g, 'depth').values())))
print 'Maximum length {0}, median depth {1}'.format(max_L, med_D)

print 'Rescaling edges'
for u, v, dat in g.edges_iter(data=True):
    
    w_uv = dat['weight']
    v_L = g.node[v]['length'] / 256.
    u_L = g.node[u]['length'] / 256.
    v_D = g.node[v]['depth'] / med_D
    u_D = g.node[u]['depth'] / med_D
    
    #beitel
    #dat['weight'] = w_uv * max_L**2 / float(u_L*v_L)

    dat['weight'] = w_uv / float(v_L*v_D * u_L*u_D)**0.5
    
nx.write_graphml(g, args.output)
