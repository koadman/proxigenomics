#!/usr/bin/env python
import gzip
import sys
import os
from Bio import SeqIO
from Bio import SeqRecord

if len(args) != 4:
    print 'Usage: [id] [input seq] [output seq]'
    sys.exit(1)

gap=SeqRecord.SeqRecord('N'*40)

with gzip.open(sys.argv[2], 'rb') as in_h:
    seq = None;
    for si in SeqIO.parse(in_h, 'fasta'):
        if not seq:
            seq = si
        else:
            seq += gap + si

seq.id=sys.argv[1]
seq.name=seq.id
seq.description='Concatentation of {0}'.format(os.path.basename(sys.argv[2]))

SeqIO.write(seq, sys.argv[3], 'fasta')
