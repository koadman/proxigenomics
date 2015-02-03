#!/usr/bin/env python
#
# Convert a SAM file to edge CSV file suitable for importing into Gephi
#
# There is an assumption that read names contain 'fwd' or 'rev' as a suffix.
#
# Eg. frg01fwd or frg9999rev
#
import sys
import math

if len(sys.argv) != 5:
    print 'Usage: [SAM file] [WGS2SCF-IDXSTATS] [edges] [nodes]'
    sys.exit(1)

class Edge:
    """Represents an edge in the network of contigs linked
    by Hi-C read pairs.
    """
    def __init__(self,nodes=[]):
        nodes.sort()
        self.nodes = nodes
        self.weight = 1
    def __eq__(self,other):
        if not isinstance(other, self.__class__):
            return False
        return self.nodes[0] == other.nodes[0] and self.nodes[1] == other.nodes[1]
    def __ne__(self,other):
        return not self.__eq__(other)
    def __str__(self):
        return '{source} {target} {weight} {normweight}'.format(source=self.nodes[0].id, target=self.nodes[1].id, weight=str(self.weight), normweight=str(self.normWeight()))
    def incWeight(self):
        self.weight += 1
    def getId(self):
        return self.nodes[0].id + self.nodes[1].id
    def normWeight(self):
        return self.weight / math.sqrt(self.nodes[0].length * self.nodes[1].length)

class Node:
    """Represents a node in the network
    """
    def __init__(self,id,length,reads):
        self.id = id
        self.length = length
        self.reads = reads
    def __eq__(self,other):
        if not isinstance(other, self.__class__):
            return False
        return self.id == other.id
    def __ne__(self,other):
        return not self.__eq__(other)
    def __lt__(self,other):
        return self.id < other.id
    def __gt__(self):
		return self.id > other.id
    def __le__(self):
		return self.id <= other.id
    def __ge__(self):
		return self.id >= other.id
    def __str__(self):
        return '{id} {length} {reads}'.format(id=self.id,length=self.length,reads=self.reads)

def updateLinkageMap(l):
    """Parse the line for new information about contig linkages. These
    may be self-self linkages or between inter-contig.
    """
    field = l.rstrip('\n').lstrip().split()
    read = field[0][:-3]
    rdir = field[0][-3:]
    contig = field[2]
    linkage = linkageMap.get(read)
    if linkage is None: linkageMap[read] = [(contig,rdir)]
    else: 
        linkage.append((contig,rdir))

# Filter lines beginning with '@' and any line where the 
# subject sequence is listed as '*'.
def filter(line):
	if line.startswith('@'): return True
	fields = line.rsplit()
	if fields[2] == '*': return True
	return False

# Read the sam file and build a linkage map
linkageMap = {}
hIn = open(sys.argv[1],'r')
[updateLinkageMap(line) for line in hIn if not filter(line)]
hIn.close()

def updateNodeMap(line):
	fields = line.rsplit()
	if len(fields) != 4:
		raise Exception('Invalid line in idxstats file [' + line + ']')
	nodeMap[fields[0]] = Node(fields[0],int(fields[1]),int(fields[2]))

# Read the idxstats file and build the node list
nodeMap = {}
hIn = open(sys.argv[2],'r')
[updateNodeMap(line) for line in hIn if not line.startswith('*')]
hIn.close()

# From the set of all linkages, convert this information
# into inter-contig edges, where the nodes are contigs.
# Count the number of redundant links as a raw edge weight.
edgeMap = {}
for (insert,linkage) in linkageMap.iteritems():
    for i in range(len(linkage)):
        for j in range(i):
            if linkage[i][1] != linkage[j][1] and linkage[i][0] != linkage[j][0]:
                e = Edge( [ nodeMap[linkage[i][0]], nodeMap[linkage[j][0]] ] )
                if e.getId() not in edgeMap:
                    edgeMap[e.getId()] = e
                else:
                    e = edgeMap.get(e.getId())
                    e.incWeight()

# Write out all edges to stdout
# Extra columns for the Gephi
hOut = open(sys.argv[3],'w')
print >>hOut, "SOURCE TARGET RAWWEIGHT WEIGHT TYPE"
for e in edgeMap.values():
    print >>hOut, e,"UNDIRECTED"
hOut.close()

hOut = open(sys.argv[4],'w')
print >>hOut, "ID LENGTH READS"
for n in nodeMap.values():
	print >>hOut, n
hOut.close()
