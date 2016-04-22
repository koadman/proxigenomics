#!/usr/bin/env python
from Bio import SeqIO
import truthtable as tt
import argparse
import sys



parser = argparse.ArgumentParser(description='Report major clusters')
parser.add_argument('-f', '--fmt', default='mcl', choices=['yaml', 'mcl'],
                    help='Input clustering format')
parser.add_argument('-m','--min', type=float, metavar='FLOAT', default=0.01,
                    help='Minimum relative size (0.01)')
parser.add_argument('--fasta',metavar='FILE',help='Fasta file of related sequences')
parser.add_argument('-v','--verbose',action='store_true',default=False,help='Verbose output')
parser.add_argument('input', metavar='FILE', help='Input clustering file')
args = parser.parse_args()

if args.fasta:
    with open(args.fasta, 'r') as fasta_h:
        seqlen = {s.id: len(s) for s in SeqIO.parse(fasta_h, 'fasta')}
    if args.verbose:
        print 'Read {0} sequence lengths'.format(len(seqlen))

if args.fmt == 'mcl':
    cl = tt.read_mcl(args.input)
elif args.fmt == 'yaml':
    cl = tt.read_yaml(args.input)
else:
    print 'Unknown input format'
    sys.exit(1)

if args.verbose:
    cl.print_tally()

n_raw = cl.num_symbols()
cl.filter_extent(args.min, seqlen)
n_filt = cl.num_symbols()

if args.verbose:
    cl.print_tally()

print n_raw, n_filt

