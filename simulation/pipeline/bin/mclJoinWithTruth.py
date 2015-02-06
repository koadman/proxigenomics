#!/usr/bin/env python

import sys

if len(sys.argv) != 4:
    print "Usage: [truth table] [clustering] [output]"
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

# Read in the clustering
with open(sys.argv[2], 'r') as h_in:
    cluster = {}
    for i, line in enumerate(h_in):
        cluster[i + 1] = line.rstrip().split()

print "Cluster table contains {0} assignments".format(len(cluster))

# Write out joined result
with open(sys.argv[3], 'w') as h_out:
    h_out.write("ctg truth predict\n")
    for cl in sorted(cluster.keys()):
        for member in sorted(cluster[cl]):
            if member in truth:
                h_out.write("{0} {1} c{2}\n".format(member, truth[member], cl))
