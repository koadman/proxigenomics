#!/usr/bin/env python
import sys
import argparse
from collections import Counter

def metisCl_to_mcl(node_map, cluster_file, output, grabbag=False):
    node_map = {}
    for line in args.node_map[0]:
        tok = line.rstrip().split()
        if len(tok) != 2:
            sys.stderr.write('Error: invalid line did not contain 2 fields [{0}]\n'.format(line))
        node_map[int(tok[0])] = tok[1]

    # create a dictionary of scaffold => cluster
    cl_set = Counter()
    isolated_nodes = []
    node_to_cl = {}
    for n, line in enumerate(args.cluster_file[0], start=1):
        tok = [int(ti) for ti in line.rstrip().split()]
        if len(tok) == 0:
            isolated_nodes.append(node_map[n])
        else:
            node_to_cl[node_map[n]] = tok
            cl_set.update(tok)

    print 'Clustering breakdown: {0}'.format(cl_set)

    # invert the dictionary, to cluster => scaffold
    cl_to_node = dict.fromkeys(cl_set)
    for k, v in node_to_cl.iteritems():
        for clid in v:
            cl = cl_to_node.get(clid)
            if cl is None:
                cl = list()
                cl_to_node[clid] = cl
            cl.append(k)

    # insert nodes which were unassigned as each being in their own cluster or in a single grab-bag
    nid = max(cl_set.keys()) + 1
    if args.grabbag:
        cl_to_node[nid] = isolated_nodes
    else:
        for node in isolated_nodes:
            cl_to_node[nid] = [node]
            nid += 1

    # sort clusters by descending size
    cl2scf = sorted(cl_to_node.values(), key=len, reverse=True)

    # write to output stream
    for i in range(0, len(cl2scf)):
        args.output.write(' '.join(cl2scf[i]) + '\n')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert Metis-style clustering file to MCL format')
    parser.add_argument('-g','--grabbag', action='store_true', default=False,
                        help='Place unassigned nodes in a single grab-bag')
    parser.add_argument('node_map', type=argparse.FileType('r'), nargs=1,
                        help='node-map file, implicit ids to names')
    parser.add_argument('cluster_file', type=argparse.FileType('r'), nargs=1,
                        help='Clustering result file')
    parser.add_argument('output', type=argparse.FileType('w'), nargs='?', default=sys.stdout,
                        help='Output file')
    args = parser.parse_args()

    metisCl_to_mcl(*vars(args))