#!/usr/bin/env python
from Bio import SeqIO
import numpy as np
import argparse
import sys


parser = argparse.ArgumentParser(description='Calculate shared bp in an N-way mauve alignment')
parser.add_argument('backbone', metavar='BACKBONE', help='Mauve backbone file')
parser.add_argument('fasta', metavar='FASTA', help='Multi-fasta used in Mauve alignment')
args = parser.parse_args()


# read fasta only to obtain sequence lengths.
seq_len = np.zeros(4, dtype=np.int32)
for i, seq in enumerate(SeqIO.parse(args.fasta, 'fasta')):
    seq_len[i] = len(seq)

# parse backbone, summing the shared columns in N, N-1, ... 1 ways

shared = None
with open(args.backbone, 'r') as bb_h:
    for line in bb_h:
        if line.startswith('seq'):
            # header line just tells us how many sequences.
            n_fields = len(line.split('\t'))
            shared = np.zeros((n_fields/2 + 1, n_fields/2))
            break

    for line in bb_h:

        line = line.strip()
        if not line:
            continue

        fields = line.split('\t')
        if len(fields) != n_fields:
            raise RuntimeError('Invalid number of fields for line [{0}], there should be {1}'.format(line, n_fields))

        fields = np.array(fields, dtype=np.int32)

        # Numpy tricks
        #   get only the non-zero columns, indicating shared extents.
        #   add their values to the shared array.
        idx = np.argwhere(fields > 0)
        n_way = len(idx) / 2
        shared[n_way, idx[::2] / 2] += np.abs(fields[idx[1::2]] - fields[idx[::2]])

    # attribute any non-reported columns to being unaligned.
    # columns of array will now sum to 1.
    shared[0, ] = seq_len - np.sum(shared[1:, ], axis=0, dtype=np.int32)

print 'Shared'
axis_sum = np.sum(shared, axis=0, dtype=np.int32)
shared = np.vstack((shared, axis_sum))
axis_sum = np.sum(shared, axis=1, dtype=np.int32)
shared = np.hstack((shared, axis_sum[:, np.newaxis]))
np.savetxt(sys.stdout, shared, fmt='%d', delimiter='\t')
print

print 'Proportion'
prop_shared = shared / np.hstack((seq_len, seq_len.sum()))
np.savetxt(sys.stdout, prop_shared, fmt='%7.5f', delimiter='\t')
print

margin = prop_shared[:, -1]
sys.stdout.write('Totals\t')
np.savetxt(sys.stdout, margin[np.newaxis,:], fmt='%7.5f', delimiter='\t')
