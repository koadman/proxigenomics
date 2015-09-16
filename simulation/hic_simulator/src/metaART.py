#!/usr/bin/env python
from Bio import SeqIO
from subprocess import call

import argparse
import sys

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Simulate a metagenomic data set from an abundance profile')
    parser.add_argument('-t', '--community-table', dest='comm_table', required=True,
                        help='Community profile table', metavar='FILE')
    parser.add_argument('-M', '--max-coverage', metavar='INT', type=int, required=True,
                        help='Coverage of must abundant taxon')
    parser.add_argument('-S', '--seed', metavar='INT', type=int, required=True, help='Random seed')
    parser.add_argument('-l', '--read-len', metavar='INT', type=int, required=True, help='Read length')
    parser.add_argument('-m', '--insert-len', metavar='INT', type=int, required=True, help='Insert length')
    parser.add_argument('-s', '--insert-sd', metavar='INT', type=int, required=True, help='Insert standard deviation')
    parser.add_argument('--art-path', default='ART_illumina', help='Path to ART executable [default: ART_illumina]')
    parser.add_argument('fasta', metavar='MULTIFASTA',
                        help='Input multi-fasta of all sequences')
    parser.add_argument('output', metavar='OUTPUT',
                        help='Output file name')
    args = parser.parse_args()
#        $ARTEXE -p -rs $SEED -m $INSERT_LEN -s $INSERT_SD -l $READ_LEN -f $X_FOLD -i $REF_SEQ -o $OUT_BASE

    profile = {}
    with open(args.comm_table, 'r') as h_table:
        for line in h_table:
            line = line.rstrip().lstrip()
            if line.startswith('#') or len(line) == 0:
                continue
            field = line.split()
            if len(field) != 3:
                print 'sequence table has missing fields at [', line, ']'
                sys.exit(1)
            profile[field[0]] = float(field[2])

    print profile

    seq_index = SeqIO.index(args.fasta, 'fasta')
    try:
        for n, seq_id in enumerate(profile):
            coverage = float(profile[seq_id] * args.max_coverage)
            print 'Generating {0} coverage for {1}'.format(coverage, seq_id)
            seq = seq_index[seq_id]
            SeqIO.write([seq], 'tmp.seq', 'fasta')
            print [args.art_path, '-p', '-rs', args.seed,
                  '-m', args.insert_len,
                  '-s', args.insert_sd,
                  '-l', args.read_len,
                  '-f', None,
                  '-i', 'tmp_seq',
                  '-o', 'reads.{0}.fq'.format(n)]
            # run ART
            # delete tmp.seq
    finally:
        if seq_index:
            seq_index.close()
