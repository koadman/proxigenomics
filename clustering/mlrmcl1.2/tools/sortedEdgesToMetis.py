
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
import os
import re

# input file has one edge per line, with each line 
# of the format "fromNode toNode edgeWeight". The edgeWeight is
# optional. The lines must be sorted by the first column.

args=sys.argv
if len(args) < 4:
	sys.stderr.write("Usage: python " + args[0] + " <inputFile>"
		+ " <numNodes> <weights?> [<numEdges>] (output on stdout)\n")
	sys.exit(-1)

regexp = re.compile("\\s+")
numNodes = args[2]
#addInLinks = int(args[3])
addInLinks = 0
weights = int(args[3])

if addInLinks > 0:
	addInLinks = True
	if weights > 0:
		sys.stderr.write("addInLinks not implemented yet" \
		+ "when weights are present.\n")
		sys.exit(-1)
else:
	addInLinks = False

if weights > 0:
	weights = True
else:
	weights = False

adjListMap = {}
wgtMap = {}

if len(args) < 6:
	from wordCount import wordcount
	numEdges=wordcount(args[1])
else:
	numEdges=int(args[5])

sys.stderr.write("Number of edges " + str(numEdges) + "\n")

#if addInLinks==True:
#	if numEdges%2 == 1:
#		sys.stderr.write("No. of edges in file " + str(numEdges) + " not even.\n")
#		sys.exit(-1)
#	else:
#		numEdges = numEdges / 2

if weights:
	print numNodes + " " + str(numEdges) + " 1"
else:
	print numNodes + " " + str(numEdges)

inputFile = file(args[1],"r")
current_row_id = 0
delimiter = " "

numPrunedEdges = 0
for line in inputFile:
	tokens = regexp.split(line.strip())
	row_id = int(float(tokens[0]))
	col_id = int(float(tokens[1]))
	if weights:
#		wt = int(tokens[2])
		wt = float(tokens[2])
#		wt = int(round(float(tokens[2])*100))
#		if wt < 1:
#			numPrunedEdges = numPrunedEdges+1
#			continue
#	print "               " + str(row_id) + "       " +	str(current_row_id) + "    " + str(col_id)

	if row_id != current_row_id:
#		print "   " + str(row_id)
		if current_row_id != 0:
			#sys.stdout.write(str(current_row_id) + " ")
			s = ' '.join(outLine)
			if weights == False and addInLinks and current_row_id in adjListMap:
				s = s + " " + ' '.join(adjListMap[current_row_id])
			print s

		outLine = [str(col_id)]
		if weights:
			outLine.append(str(wt))
#		print "outLine:" + ' '.join(outLine)
		current_row_id = current_row_id + 1
		if current_row_id % 5000 == 0:
			sys.stderr.write("Done with " + str(current_row_id) + 
			" nodes.\n")
	else:	
		outLine.append(str(col_id))
		if weights:
			outLine.append(str(wt))

	while row_id > current_row_id:
#		print "" + str(row_id) + "  " + str(current_row_id)
		if weights==False and addInLinks and current_row_id in adjListMap:
			l = adjListMap[current_row_id]
			print ' '.join(l)
		else:
			print
		current_row_id = current_row_id + 1
		if current_row_id % 5000 == 0:
			sys.stderr.write("Done with " + str(current_row_id) + 
			" nodes.\n")

	if weights == False and addInLinks:
		if row_id not in adjListMap:
			list = [str(col_id)]
			adjListMap[row_id] = list
		else:
			adjListMap[row_id].append(str(col_id))

#sys.stdout.write( str(row_id) + " " )

s = ' '.join(outLine)
if addInLinks and row_id in adjListMap:
	s = s + " " + ' '.join(adjListMap[row_id])
print s

numNodes = int(numNodes)
while row_id < numNodes:
	row_id = row_id + 1
#	print row_id
#	print numNodes
#	print row_id < numNodes
	if addInLinks and row_id in adjListMap:
		print ' '.join(adjListMap[row_id])
	else:
		print 

