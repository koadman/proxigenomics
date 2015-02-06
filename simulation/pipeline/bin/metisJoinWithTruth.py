#!/usr/bin/env python

import sys

if len(sys.argv) != 5:
    print "Usage: [truth table] [node map] [clustering] [output]"
    sys.exit(1)

# Read in the truth table
with open(sys.argv[1], 'r') as h_in:
    truth = {}
    for line in h_in:
        tok = line.split()
        if len(tok) != 3:
            print "Invalid line in truth table, expecting three fields [{0}]".format(line)
            sys.exit(1)
        truth[tok[0]] = tok[1]

print "Truth table contains {0} assignments".format(len(set(truth.values())))

# Read node # to name mapping file
with open(sys.argv[2], 'r') as h_in:
    node = {}
    for line in h_in:
        tok = line.rstrip().split()
        if len(tok) != 2:
            raise IOError('Invalid number of fields on line {0}'.format(line))
        node[int(tok[0])] = tok[1]

# Read in the clustering
with open(sys.argv[3], 'r') as h_in:
    cluster = {}
    for i, line in enumerate(h_in):
        clid = int(line.rstrip()) + 1
        cl = cluster.get(clid)
        if cl is None:
            cl = []
            cluster[clid] = cl
        cl.append(node[i + 1])

print "Cluster table contains {0} assignments".format(len(cluster))

# Write out joined result
with open(sys.argv[4], 'w') as h_out:
    h_out.write("ctg truth predict\n")
    for cl in sorted(cluster.keys()):
        for member in sorted(cluster[cl]):
            if member in truth:
                h_out.write("{0} {1} c{2}\n".format(member, truth[member], cl))
