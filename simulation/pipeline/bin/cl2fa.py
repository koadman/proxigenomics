#!/usr/bin/env python
from Bio import SeqIO
import truthtable as tt
import os
import sys

if len(sys.argv) != 4:
    print 'usage <clustering> <multi-fassta> <out dir>'
    sys.exit(1)

seqidx = SeqIO.index(sys.argv[2], 'fasta')

lengths = {}
for si in seqidx:
    lengths[si] = len(seqidx[si])

print 'Calculated {0} sequence lengths'.format(len(lengths))

t = tt.read_mcl(sys.argv[1])
t.print_tally()

#t.filter_extent(0.01, lengths)
cl2seq = t.invert()

for clid, seqids in cl2seq.iteritems():
    seqs = []
    for si in seqids:
        seqs.append(seqidx[si])
    print 'For cluster {0} writing {1} sequences'.format(clid,len(seqs))
    SeqIO.write(seqs, os.path.join(sys.argv[3], 'cl{0}.fasta'.format(clid)), 'fasta')

