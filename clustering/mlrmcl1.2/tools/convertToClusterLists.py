 # MLR-MCL (Multi-Level Regularized Markov Clustering) - Version 1.2
 # 
 # Copyright 2010, The Ohio State University. 
 # Portions Copyright 1997, Regents of the University of Minnesota.
 # All rights reserved.
 # 
 # Redistribution and use in source and binary forms, with or
 # without modification, are
 # permitted provided that the following conditions are met:
 # 
 # 1. Redistributions, with or without modifications, of source 
 # code must retain the above
 # copyright notice, this list of conditions and the following
 # disclaimer.
 # 
 # 2. Redistributions, with or without modifications, in binary form 
 # must reproduce the above
 # copyright notice, this list of conditions and the following
 # disclaimer
 # in the documentation andor other materials provided with the
 # distribution.
 # 
 # 3. The names of the Ohio
 # State University, the University of Minnesota and
 # their contributors may not be used to endorse or promote products
 # derived from this software without specific prior permission. 
 # 
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 # CONTRIBUTORS
 # "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 # LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 # FOR
 # A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 # COPYRIGHT
 # OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 # INCIDENTAL,
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 # USE,
 # DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 # ANY
 # THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 # TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 # USE
 # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 # DAMAGE.
 # 
 #

 # Contributors: Venu Satuluri, Srinivasan Parthasarathy
 # 
 # Reference: "Scalable Graph Clustering using Stochastic Flows",
 # KDD 2009. http:doi.acm.org10.11451557019.1557101
 #


import sys
import re

# The output from MLR-MCL (as also Metis and Graclus) 
# has as many lines as nodes in the graph, and for 
# each line has the cluster no. to which the node represented by
# this line no. belongs to. This script takes such a file and
# outputs a file with one cluster per line and all the members of
# the cluster given in that line in a space separated format.

def converttoclusterlists(fname):
	clusters = {}
	inputFile = file(fname,"r")

	lineNo = 1
	reDoubleColumn = re.compile("^([0-9]+)\\s+([0-9]+)$")

	for line in inputFile:
		line = line.strip()
		if reDoubleColumn.match(line):
			m = reDoubleColumn.match(line)
			clusterNo = int(float(m.group(2)))
			nodeId = int(float(m.group(1)))
		else:
		#	clusterNo=int(line)
			clusterNo=int(float(line))
			nodeId = lineNo

		if clusterNo not in clusters:
			clusters[clusterNo]=[]
		clusters[clusterNo].append(nodeId)
		lineNo=lineNo+1

	return clusters

# start of main
if __name__ == "__main__":
	args = sys.argv
	if len(args) < 2:
		sys.stderr.write("Usage: python " + args[0] + " <inputFile> "
			+ "(output on stdout)\n")
		sys.exit(-1)

	clusters = convertmetisoutput(args[1])
	for i in clusters:
		print ' '.join(map(str,clusters[i]))

