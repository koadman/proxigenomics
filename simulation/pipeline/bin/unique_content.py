#!/usr/bin/env python
from collections import Counter
import networkx as nx
import community as com
import numpy as np
import sys

def ordered_counts(counter):
    keys = sorted(counter.keys())
    return dict((k, counter[k]) for k in keys)

def relative_proportion(counter):
    N = float(np.sum(counter.values()))
    keys = sorted(counter.keys())
    return dict((k, counter[k]/N) for k in keys)

def linkage_test(g, u, v):
    if g.has_edge(u,v):
        return 'Direct'
    try:
        pl = len(nx.shortest_path(g,v,u))-1
        return 'Indirect'
    except nx.NetworkXNoPath:
        return 'Disconnected'

if len(sys.argv) != 6:
    print 'Usage: [MIN ID] [MIN COV] [MIN LENGTH] [CTG2REF PSL] [HIC2CTG graphml]'
    sys.exit(1)

min_id = float(sys.argv[1])
min_cov = float(sys.argv[2])
min_length = int(sys.argv[3])

ctg2ref_map = {}

# Parse PSL file for uniquely mapping contigs
with open(sys.argv[4], 'r') as psl_file:
    #outh = open('perid.csv','w')
    for line in psl_file:
        line = line.strip()
        if not line:
            break
        fields = line.split()
        alen = int(fields[12]) - int(fields[11])
        qlen = int(fields[10])
        qname = fields[9]
        sname = fields[13]


        matches = int(fields[0])
        mismatches = int(fields[1])
        repmatches = int(fields[2])
        q_num_insert = int(fields[4])
        perid = (1.0 - float(mismatches + q_num_insert) / float(matches + mismatches + repmatches))
        cov = float(alen)/float(qlen)

        #outh.write('{0} {1} {2} {3}\n'.format(alen, qlen, cov, perid))

        if perid >= min_id and cov >= min_cov and qlen >= min_length:
            # get read depth from SPADES contig name    
            #toks = qname.spli8t('_')
            #rd = float(toks[5])
            #if rd < 60.0 and
            #print qname, sname, alen, qlen, cov
            if qname not in ctg2ref_map:
                ctg2ref_map[qname] = {'count': 1, 'data':{sname: (alen, qlen, cov, perid)}}
            else:
                ctg2ref_map[qname]['count'] += 1
                ctg2ref_map[qname]['data'][sname] = (alen, qlen, cov, perid)

# prune all non-unique
for k in ctg2ref_map.keys():
    if ctg2ref_map[k]['count'] != 1:
        del ctg2ref_map[k]

print 'There are {0} uniquely mapping contigs'.format(len(ctg2ref_map))
#for k, v in ctg2ref_map.iteritems():
#    print k, v

# Hi-C contig graph
print 'Reading graphml file'
g = nx.read_graphml(sys.argv[5])

print 'Random sampling'
rnd_count = Counter()
for n in xrange(5):
    print 'round {0}'.format(n)
    for ri in g.nodes():
        while True:
            rj = np.random.choice(g.nodes())
            if rj != ri:
                rnd_count[linkage_test(g, ri, rj)] += 1
                break

print 'Random', ordered_counts(rnd_count), relative_proportion(rnd_count)

unmatched = 0
counter = Counter()
unique_ids = ctg2ref_map.keys()
N_ids = len(unique_ids)
for i in xrange(N_ids):

    if N_ids > 20 and i % (N_ids /20) == 0:
        print '{0}/{1}'.format(i, N_ids)
    
    v = unique_ids[i]
    # skip if not unique or node doesn't exists in graph
    if ctg2ref_map[v]['count'] != 1 or not g.has_node(v):
        continue

    v_src = ctg2ref_map[v]['data'].keys()[0]

    for j in xrange(i+1, N_ids):

        u = unique_ids[j]
        # skip if not unique or node doesn't exists in graph
        if ctg2ref_map[u]['count'] != 1 or not g.has_node(u):
            continue

        u_src = ctg2ref_map[u]['data'].keys()[0]

        # u,v directly linked
        if g.has_edge(v, u):
            counter['Direct'] += 1
            #print u, v
            if v_src != u_src:
                unmatched += 1
                #print 'Direct with diff src: {0},{3} {1},{4} w:{2}\n{5}\n{6}'.format(v_src,u_src, g[v][u], v, u, ctg2ref_map[v], ctg2ref_map[u])
        else:
            #counter['Other'] += 1
            # find shortest path between u,v
            if nx.has_path(g, v, u):
                counter['Indirect'] += 1
            else:
                counter['Disconnected'] += 1

print 'Unique', ordered_counts(counter), relative_proportion(counter)
print 'Unmatched sources for uniques',unmatched

