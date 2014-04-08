#!/usr/bin/env python

import sys

if len(sys.argv) != 5:
    print "Usage: [truth table] [node map] [clustering] [output]"
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

# Read node # to name mapping file
hIn = open(sys.argv[2], 'r')
node = {}
for line in hIn:
	tok = line.rstrip().split()
	if len(tok) != 2:
		raise IOError('Invalid number of fields on line {0}'.format(line))
	node[int(tok[0])] = tok[1]

# Read in the clustering
hIn = open(sys.argv[3],'r')
cluster = {}
for i,line in enumerate(hIn):
	clid = int(line.rstrip())+1
	cl = cluster.get(clid)
	if cl is None:
		cl = []
		cluster[clid] = cl
	cl.append(node[i+1])
hIn.close()

print "Cluster table contains {0} assignments".format(len(cluster))

# Write out joined result
hOut = open(sys.argv[4],'w')
hOut.write("ctg truth predict\n")
for cl in sorted(cluster.keys()):
    for member in sorted(cluster[cl]):
        if member in truth:
            hOut.write("{0} {1} c{2}\n".format(member,truth[member],cl))
hOut.close()
