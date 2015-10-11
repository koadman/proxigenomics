#!/bin/env python

from Bio import SeqIO
import argparse
import os

parser = argparse.ArgumentParser(description='Split multifasta file')
parser.add_argument('-f', '--force', default=False, action='store_true', help='Overwrite existing output files')
parser.add_argument('fasta', metavar='FASTA_FILE', type=argparse.FileType('r'), help='Input multifasta file')
parser.add_argument('output_path', metavar='PATH', help='Output path')

allow_chars = ('_', '-', '.')

args = parser.parse_args()

if not os.path.exists(args.output_path):
    raise IOError('Output path {0} does not exist'.format(args.output_path))

for n, seq in enumerate(SeqIO.parse(args.fasta, 'fasta'), start=1):
    oname = '{0}.fasta'.format(''.join(ch for ch in seq.id if ch.isalnum() or ch in allow_chars).rstrip())
    oname = os.path.join(args.output_path, oname)
    if not args.force and os.path.exists(oname):
        raise IOError('An output file already exists with the name {0}'.format(oname))
    SeqIO.write(seq, oname, 'fasta')

print 'Multi-fasta {0} contained {1} sequences'.format(args.fasta.name, n)
