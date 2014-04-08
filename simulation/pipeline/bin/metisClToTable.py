#!/usr/bin/env python
import sys

if len(sys.argv) != 4:
	print 'Usage: [node map] [metis cluster] [table out]'
	sys.exit(1)

try:
	hmap = open(sys.argv[1],'r')
	hcl = open(sys.argv[2],'r')
	hout = open(sys.argv[3],'w')
	
	try:
		node_map = {}
		for line in hmap:
			tok = line.rstrip().split()
			if len(tok) != 2:
				print "Error: invalid line did not contain 2 fields [{0}]".format(line)
				raise IOError
			node_map[int(tok[0])] = tok[1]
	finally:
		hmap.close()

	try:
		cl_set = set()
		cl_assign = {}
		for i,line in enumerate(hcl):
			tok = line.rstrip().split()
			cl_assign[node_map[i+1]] = tok
			cl_set.update(tok)

		# set up a lookup table of cluster to scaffold
		cl_dict = dict.fromkeys(cl_set)
		for k,v in cl_assign.iteritems():
			for clid in v:
				cl = cl_dict.get(clid)
				if cl is None:
					cl = list()
					cl_dict[clid] = cl
				cl.append(k)
	finally:
		hcl.close()

	try:
		# sort the lookup table longest to shortest and write out
		cl2scf = sorted(cl_dict.values(), key=len, reverse=True)
		
		for i in range(0,len(cl2scf)):
			hout.write(' '.join(cl2scf[i]) + '\n')

	finally:
		hout.close()
		
except IOError, e:
	print "Error", e
