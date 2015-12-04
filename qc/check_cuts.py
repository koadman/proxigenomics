#!/usr/bin/env python
from __future__ import division
import pysam
import sys

fname = sys.argv[1]
cut_site = sys.argv[2]
track = pysam.AlignmentFile(fname, "rb")
count = 0
found_count = 0
for aln in track.fetch(until_eof=True):
    if(aln.mapping_quality < 20):
        continue    
    # skip alignments that are too close to the contig boundaries
    if(aln.reference_start <= 5 or aln.reference_start + 150 >= track.header['SQ'][aln.reference_id]['LN']):
        continue
    if(len(aln.cigar) > 1):
        pos = 0
        for chunk in aln.cigar:
            if(chunk[0] == 0):
                if(pos > 4):
                    print aln.query_name
                    print aln.cigar
                    print "from " + str(pos-4)
                    ss = aln.query_sequence[ (pos-4):(pos+6) ]
                    print ss
                    if(ss.find(cut_site) >= 0):
                        found_count += 1
                    count += 1
                if(pos + chunk[1] < len(aln.query_sequence) - 4):
                    print aln.query_name
                    print aln.cigar
                    print "from2 " + str(pos + chunk[1]-6)
                    ss = aln.query_sequence[ (pos + chunk[1]-6):(pos + chunk[1] + 2) ]
                    print ss
                    if(ss.find(cut_site) >= 0):
                        found_count += 1
                    count += 1
            pos = pos + chunk[1]

print str(found_count) + "/" + str(count) + " = " + str(found_count / count)

