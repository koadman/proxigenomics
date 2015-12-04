#!/usr/bin/env python
# remove sequences with long mapping regions
from __future__ import division
import pysam
import sys

fname = sys.argv[1]
track = pysam.AlignmentFile(fname, "rb")
prevmapped = -1
prev_aln = None
for aln in track.fetch(until_eof=True):
    if(prevmapped == -1):
        prev_aln = aln
    mapped = 0
    total = 0
    for chunk in aln.cigar:
        if(chunk[0] == 0):
            total += chunk[1]
    if(total > 40):
        mapped = 1
    if(prev_aln.query_name == aln.query_name and mapped == 0 and prevmapped == 0):
        print "@" + prev_aln.query_name
        print prev_aln.query_sequence
        print "+"
        print pysam.toQualityString(prev_aln.query_qualities)
        print "@" + aln.query_name
        print aln.query_sequence
        print "+"
        print pysam.toQualityString(aln.query_qualities)

    prev_aln = aln
    prevmapped = mapped


