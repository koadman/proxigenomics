
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
from convertToClusterLists import converttoclusterlists

permArrays = {}
factArray = [0, 0]

def log10fact(n):
	if n < 0:
		raise ValueError('cannot calculate factorial of n for n='+str(n)+'!')
	n = int(n)
	if len(factArray) <= n:
		for i in range(len(factArray),n+1):
			factArray.append(factArray[i-1] + math.log10(float(i)))
	
	return factArray[n]
	

def log10Perm(n,x):
	if x > n:
		raise ValueError('cannot calculate n perm x for n='+str(n)+', x='+str(x))

	intn = int(n)
	intnx = int(n-x)
	intx = int(x)
	if intx == 1:
		return math.log10(float(n))
	if intn not in permArrays:
		permArrays[intn] = []
		permArrays[intn].append(float(0))
	if intx not in permArrays[intn]:
		for i in range(len(permArrays[intn]),intx+1):
			permArrays[intn].append(permArrays[intn][i-1] +	math.log10(float(n-i+1)))
	return permArrays[intn][intx]

#	result = float(0)
#	for i in range(intn,intnx,-1):
#		result = result + math.log10(float(i))
#
#	return result


# log base 10 of (n choose x) 
def log10BinCoeff(n,x):
	if x > n or x < 0:
		raise ValueError('cannot calculate n choose x for n='+str(n)+', x='+str(x))

	if x > n/2:
		x = n-x;
	
	#result = log10Perm(n,x)-log10Perm(x,x)
	result = log10Perm(n,x)-log10fact(x)
	return result


args = sys.argv

if len(args) < 7:
	sys.stderr.write("Usage: python " + args[0] + \
	" <outputClusteringFile> <groundTruthClusteringFile> " + \
	"<groundTruthSizesFile> <minSize> <noOfElements> " + \
	"<maxSize> [<fscores>?] [<normalizer>] [<intersectAtLeast2>] " + \
	"(output on stdout)\n")
	sys.exit(-1)

#outputClustering file can either single column or double column.
#If 2 columns are present, it is assumed that first is the element
#id and the second is the cluster id; in the case of 1 column, the
#line no is assumed to be the element id.
#
#groundTruthClustering must be double
#column, with similar semantics as above. The file should be
#sorted by the first column. Each element can belong
#to multiple clusters.
#
#groundTruthSizes is double column, with first column being size
#of the cluster followed by the cluster id.
#
#minSize is the minimum size a cluster needs to have before we
#evaluate it (or use it in the evaluation, i.e. it applies to both
#output clusters as well as ground truth clusters).

#outClusters has output clusters i.e. clusters to be evaluated
#it's a map from clusterId to the list of elements in that
#clusterId.
outClusters = converttoclusterlists(args[1])
sys.stderr.write("Finished loading file " + args[1] + "\n")

# sizes of gt clusters, indexed by the cluster id.
gtClusterSizes = {}

# ground truth memberships, stored as a list of lists, one list
# per node id. Nodes that don't belong to any cluster are given
# an empty list.
gtMemberships = []

minSize = int(args[4])
maxSize = int(args[6])
N = float(args[5])
isFvalues = False 
inputNormalizer = 0
intersect2 = False
if len(args) > 7:
	if args[7].startswith("y") or args[7].startswith("Y") or args[7].startswith("1"):
		isFvalues = True
if len(args) > 8:
	inputNormalizer = float(args[8])
if len(args) > 9:
	if args[9].startswith("y") or args[9].startswith("Y") or args[9].startswith("1"):
		intersect2 = True


rx=re.compile("\\s+")

# first populate gtClusterSizes
for line in file(args[3], "r"):
	line = line.strip()
	tokens = rx.split(line)
	size = int(tokens[0])
	if size >= minSize and size <= maxSize:
		gtClusterSizes[int(tokens[1])] = size

sys.stderr.write("Finished loading file " + args[3] + "\n")
sys.stderr.write(str(len(gtClusterSizes)) + \
" ground truth clusters with size >= " + str(minSize) + \
" and <= " + str(maxSize) + "\n")

# now populate gtMemberships
for line in file(args[2], "r"):
	line = line.strip()
	tokens = rx.split(line)
	nodeId = int(tokens[0])
	while len(gtMemberships) < nodeId:
		gtMemberships.append([])

	clusterId = int(tokens[1])
	if clusterId in gtClusterSizes:
		gtMemberships[nodeId-1].append(clusterId)

while len(gtMemberships) < int(N):
	gtMemberships.append([])

#for i in gtMemberships:
#	print ' '.join(map(str,i))

sys.stderr.write("Finished loading file " + args[2] + "\n")

wgtdAvgF = float(0)
wgtdAvgP = float(0)
wgtdAvgR = float(0)
sumPvalues = float(0)
sumSquaresPvalues = float(0)

bestPvalue = 0
pvalueThreshold = 2 #significance threshold is 0.01

NbinCoeffs = {} 
#a map for storing (N choose x) coefficients, to
#avoid expensive recomputation

# now go through each output cluster and compute the best f-value
# for each cluster.
sumOfMs = 0
actualNumOutClusters = 0
print '#Scroll to bottom for final scores'
print '#predictedClusterId groundTruthClusterId Num_Intersect F-score Precision Recall'
#print '%d %d %d %.3f %.3f %.3f'%(clusterId,	bestf_gtClusterId, M, bestf, bestp, bestr)
for clusterId in outClusters:
	cluster = outClusters[clusterId]
	sumOfMs = sumOfMs + len(cluster)
#	print "---------"
#	print ' '.join(map(str,cluster))
	if len(cluster) < minSize or len(cluster) > maxSize:
		continue
	else:
		actualNumOutClusters = actualNumOutClusters + 1
	
	# store the size of intersection of this cluster with each
	# ground truth cluster it has a non-zero intersection with.
	intersects = {}
	for i in cluster:
		# i-1 'cause gtMemberships is a list and not an
		# associative array.
			
		for gtCluster in gtMemberships[i-1]:
			if gtCluster not in intersects:
				intersects[gtCluster] = 1
			else:
				intersects[gtCluster] = intersects[gtCluster] + 1
#			if clusterId == 1:
#				print str(gtCluster)+" "+str(intersects[gtCluster])
	
#	print "------------"
#	for k in intersects:
#		print str(k)+" "+str(intersects[k])

#	print len(intersects)

	#done calculating size of intersections
	#calculate best f-score
	bestf = float(0)
	bestp = float(0)
	bestr = float(0)
	bestf_gtClusterId = -1
	bestpvalue = float(0)
	bestpvalue_gtClusterId = -1

	M = float(len(cluster))
	for gtCluster in intersects:
		n = float(gtClusterSizes[gtCluster])
		m = float(intersects[gtCluster])
		prec = m/M
		rec = m/n
		fmeasure = 2*prec*rec/(prec+rec)
		if m > M or m < 0:
			sys.stderr.write("Yikes! m (" + str(m) + ") > M"
			+ "(" + str(M) + ")")
			sys.exit(-1)
		if isFvalues:
			if fmeasure > bestf:
				if intersect2 == False or \
				(intersect2 == True and m > 1):
					bestf = fmeasure
					bestf_gtClusterId = gtCluster
					bestp = prec
					bestr = rec
		else:
			if M in NbinCoeffs:
				pvalue = NbinCoeffs[M]
			elif N-M in NbinCoeffs:
				pvalue = NbinCoeffs[N-M]
			else:
				pvalue = log10BinCoeff(N,M)
				NbinCoeffs[M] = pvalue
			subtract = log10BinCoeff(n,m) +	log10BinCoeff(N-n,M-m)	
# we are trying to approximate the sum of probabilities with the
# biggest probability. 
#			for i in range(int(m),min(int(n),int(M))+1):
#				t = log10BinCoeff(n,i) + log10BinCoeff(N-n,M-i)
#				if t > subtract:
#					subtract = t
			pvalue = pvalue - subtract
#			pvalue = pvalue - log10BinCoeff(n,m)
#			pvalue = pvalue - log10BinCoeff(N-n,M-m)
			if pvalue > bestpvalue:
				bestpvalue = pvalue
				bestpvalue_gtClusterId = gtCluster
				bestf = fmeasure
				bestp = prec
				bestr = rec

	sumPvalues = sumPvalues + bestpvalue
	sumSquaresPvalues = sumSquaresPvalues +	bestpvalue*bestpvalue
	wgtdAvgF = wgtdAvgF + M*bestf
	wgtdAvgP = wgtdAvgP + M*bestp
	wgtdAvgR = wgtdAvgR + M*bestr
	if isFvalues and bestf > 0:
		print '%d %d %d %.3f %.3f %.3f'%(clusterId,	bestf_gtClusterId, M, bestf, bestp, bestr)
#		print '{0:d} {2:d} {1:.3f} {3:.3f} {4:.3f}'.format(clusterId, bestf_gtClusterId, bestf, bestp, bestr)
#		print(str(clusterId)+" "+str(bestf)+" "+str(bestf_gtClusterId))
	elif bestpvalue > pvalueThreshold:
		print '%d %d %.3f %.3f %.3f %.3f'%(clusterId,bestpvalue_gtClusterId, bestpvalue, bestf, bestp, bestr)

print 'Num_Clusters: %d'%(len(outClusters))
print 'Num_Clusters >= %d and <= %d: %d'%(minSize, maxSize,
actualNumOutClusters)

#if isFvalues:
#for normalizer in [sumOfMs, N]:
#for normalizer in [sumOfMs]:
if inputNormalizer > 0:
	norm = inputNormalizer
#	norm = sumOfMs
else:
	norm = sumOfMs
for normalizer in [norm]:
	wavgF = wgtdAvgF / normalizer
	wavgP = wgtdAvgP / normalizer
	wavgR = wgtdAvgR / normalizer
#	print 'Normalizing by %d'%(normalizer)
	print 'Avg F-score: %.4f'%(wavgF)
	print 'Avg Precision: %.4f'%(wavgP)
	print 'Avg Recall: %.4f'%(wavgR)
	print '(Note the Avg F-score is not necessarily the harmonic mean of the Avg Precision and Avg Recall)'


if isFvalues == False:
	l = float(len(outClusters))
	meanp = sumPvalues/l
	print 'Mean P-value:%.4f'%(meanp)	
	stdDev = l*sumSquaresPvalues - sumPvalues*sumPvalues
	stdDev = stdDev / (l*(l-1))
	stdDev = math.sqrt(stdDev)
	print 'Std. Devn. of P-values:%.4f'%(stdDev)

