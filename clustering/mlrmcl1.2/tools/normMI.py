
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
import math

# Calculate Normalized Mutual Information (NMI)
# This program accepts two files as input. First is the file with
# the ground truth assignments of ids to clusters (classes). The
# first column should be id and the second should be class. The
# second file contains predicted cluster assignments. Again,
# first column is id and second is cluster assignments.
# NOTE: This program may not give the same results when the order
# of the input files is changed and the ground truth is not
# present for all the ids that are assigned clusters in the
# clusterAssignment file. This is because the size of the ground
# truth cluster is different if the cluster assignment file is
# processed first. 

args = sys.argv
if len(args) < 3:
	sys.stderr.write("Usage: python " + args[0] + " <groundTruthClusterFile> "
		+ " <predictedClusterFile> (output on stdout)\n")
	sys.exit(-1)

clusterMap = {}
gtMap = {}
gtCounts = []
clusterCounts = []
gtClusterCounts = {}
gtAssignments = {}
regexp = re.compile("\\s+")

lineNo = 1
for line in file(args[1],"r"):
	t = regexp.split(line.strip())
	if len(t) > 1:
		id = int(t[0])
		gtc = int(t[1])
	else:
		gtc = int(t[0])
		id = lineNo

	if gtc not in gtMap:
		gtMap[gtc] = len(gtMap)
		gtCounts.append(0)
	gtCounts[gtMap[gtc]] = gtCounts[gtMap[gtc]]+1
	gtAssignments[id] = gtMap[gtc]
	lineNo = lineNo+1

lineNo = 1
N_cluster = 0
for line in file(args[2],"r"):
	t = regexp.split(line.strip())
	if len(t) > 1:
		id = int(t[0])
		cass = int(t[1])
	else:
		id = lineNo
		cass = int(t[0])

	if id in gtAssignments:
		if cass not in clusterMap:
			clusterMap[cass] = len(clusterMap)
			clusterCounts.append(0)
		cmapass = clusterMap[cass]
		clusterCounts[cmapass] = clusterCounts[cmapass]+1
		N_cluster = N_cluster+1

		gtass = gtAssignments[id]
		if gtass not in gtClusterCounts:
			gtClusterCounts[gtass] = {}
		if cmapass not in gtClusterCounts[gtass]:
			gtClusterCounts[gtass][cmapass] = 0
		gtClusterCounts[gtass][cmapass] = gtClusterCounts[gtass][cmapass]+1

	lineNo = lineNo+1

mi = 0 # mutual information
N_gt = float(len(gtAssignments))
N = N_gt
#print N
#print N_cluster
logN = math.log(N,2)
logN_cluster = math.log(float(N_cluster),2)
h_gtc = 0 # entropy of ground truth assignments
h_cluster = 0 #entropy of cluster assignments

for gtc in gtClusterCounts:
	gtcc = gtClusterCounts[gtc]
	gtcCount = math.log(float(gtCounts[gtc]),2)
	for cluster in gtcc:
		clusterCount = float(clusterCounts[cluster])
		gtccCount = float(gtcc[cluster])
		t = gtccCount/N
		u = logN+math.log(gtccCount,2)-gtcCount-math.log(clusterCount,2)
		mi = mi+t*u

for cCount in clusterCounts:
	cCount = float(cCount)
	h_cluster = h_cluster+(cCount/N_cluster)*(logN_cluster-math.log(cCount,2))

for gtCount in gtCounts:
	gtCount = float(gtCount)
	h_gtc = h_gtc+(gtCount/N)*(logN-math.log(gtCount,2))

nmi = mi*2/(h_gtc+h_cluster)
print 'Norm. Mutual Information: %.4f'%(nmi)
print 'Mutual Information: %.4f'%(mi)
print 'Entropy of Ground Truth Clustering: %.4f'%(h_gtc)
print 'Entropy of Predicted Clustering: %.4f'%(h_cluster)
#print str(nmi) + ' ' + str(mi) + ' ' + str(h_gtc) + ' ' + str(h_cluster)
		
