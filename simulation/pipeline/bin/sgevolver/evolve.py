#!/usr/bin/env python

from Bio import SeqIO, Phylo
import subprocess
import argparse
import os.path

parser = argparse.ArgumentParser(description='Evolve a set of sequences from a starting reference')
parser.add_argument('--anc-start', metavar='INT', default=0, type=int, help='Ancestral sequence start on input')
parser.add_argument('--work', metavar='PATH', default='.', help='Working directory path')
parser.add_argument('length', type=int, metavar='LENGTH', nargs=1, help='Evolved sequence length (bp)')
parser.add_argument('tree', metavar='NEWICK', nargs=1, help='Evolutionary tree in Newick format')
parser.add_argument('sequence', metavar='FASTA', nargs=1, help='Initial sequence in FASTA format')
args = parser.parse_args()

with open(args.sequence[0], 'r') as h_in:
    input_seq = SeqIO.read(h_in, 'fasta')

with open(args.tree[0], 'r') as h_in:
    tree = Phylo.read(h_in, 'newick')

start = args.anc_start
stop = start + args.length[0]
if stop > len(input_seq):
    raise Exception('Resulting endpoint of ancestral sequence goes beyond its maximum length {0}'
                    .format(len(input_seq)))

anc_length = stop - start
with open(os.path.join(args.work, 'ancestral_donor.raw'), 'w') as h_out:
    print 'Writing ancestral sequence of {0}bp beginning at {1}'.format(anc_length, start)
    h_out.write(str(input_seq[start:stop].seq))
    h_out.write('\n')


