#!/usr/bin/env python

import sys
import yaml

if len(sys.argv) != 4:
    print "Usage: [truth table] [clustering] [output]"
    sys.exit(1)

# Read in the truth table
truth = None
with open(sys.argv[1], 'r') as h_in:
    truth = yaml.load(h_in)

print "Truth table contains {0} objects assigned to {1} classes".format(len(truth), len(set(truth.values())))

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
