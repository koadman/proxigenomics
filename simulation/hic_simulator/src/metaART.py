#!/usr/bin/env python
from Bio import SeqIO

import argparse
import os
import subprocess
import sys

TMP_INPUT = 'seq.tmp'
TMP_OUTPUT = 'reads.tmp'

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
    parser.add_argument('output_base', metavar='OUTPUT BASE',
                        help='Output file name')
    args = parser.parse_args()

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
    with open('{0}1.fq'.format(args.output_base), 'w') as output_R1, \
            open('{0}2.fq'.format(args.output_base), 'w') as output_R2:
        try:
            for seq_id in profile:
                coverage = float(profile[seq_id] * args.max_coverage)
                print 'Generating {0} coverage for {1}'.format(coverage, seq_id)
                seq = seq_index[seq_id]
                SeqIO.write([seq], TMP_INPUT, 'fasta')
                subprocess.call([args.art_path,
                                 '-p',   # paired-end sequencing
                                 '-na',  # no alignment file
                                 '-rs', str(args.seed),
                                 '-m', str(args.insert_len),
                                 '-s', str(args.insert_sd),
                                 '-l', str(args.read_len),
                                 '-f', str(coverage),
                                 '-i', TMP_INPUT,
                                 '-o', TMP_OUTPUT])
                with open('{0}1.fq'.format(TMP_OUTPUT), 'r') as tmp_h:
                    output_R1.write(tmp_h.read())
                    os.remove(tmp_h.name)
                with open('{0}2.fq'.format(TMP_OUTPUT), 'r') as tmp_h:
                    output_R2.write(tmp_h.read())
                    os.remove(tmp_h.name)
                os.remove(TMP_INPUT)
        finally:
            if seq_index:
                seq_index.close()
