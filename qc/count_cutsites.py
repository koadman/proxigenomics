#!/usr/bin/env python
from __future__ import division
import sys
import gzip

fq = gzip.open(sys.argv[1])
cutseqs = sys.argv[2:6]
print "Checking cutseqs " + " ".join(cutseqs)
i = -1
cutcount = [0,0,0,0]
totalseq = 0
for line in fq:
    i += 1
    if(i % 4 == 1):
        for j in range(len(cutseqs)):
            cutcount[j] += line.count(cutseqs[j])
    totalseq += len(line)
    if(totalseq > 100000000):
        break

denom = cutcount[0]
for j in range(len(cutseqs)):
    cutcount[j] /= denom

print " ".join(map(str,cutcount)) + " total sites: " + str(totalseq) 

