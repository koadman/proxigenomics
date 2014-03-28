#!/usr/bin/env python

import sys

if len(sys.argv) != 4:
    print "Usage: [truth table] [clustering] [output]"
    sys.exit(1)

# Read in the truth table
hIn = open(sys.argv[1],'r')
truth = {}
for line in hIn:
    tok = line.split()
    if len(tok) != 3:
        print "Invalid line in truth table, expecting three fields [{0}]".format(line)
        sys.exit(1)
    truth[tok[0]]=tok[1]
hIn.close()

print "Truth table contains {0} assignments".format(len(set(truth.values())))

# Read in the clustering
hIn = open(sys.argv[2],'r')
cluster = {}
for i,line in enumerate(hIn):
    cluster[i+1] = line.split()
hIn.close()

print "Cluster table contains {0} assignments".format(len(cluster))

# Write out joined result
hOut = open(sys.argv[3],'w')
hOut.write("ctg truth predict\n")
for cl in sorted(cluster.keys()):
    for member in sorted(cluster[cl]):
        if member in truth:
            hOut.write("{0} {1} c{2}\n".format(member,truth[member],cl))
hOut.close()
