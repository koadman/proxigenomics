#!/usr/bin/env python
from Bio import SeqIO

import argparse
import os
import subprocess
import sys

TMP_INPUT = 'seq.tmp'
TMP_OUTPUT = 'reads.tmp'
R1_FILE = '{0}1.fq'.format(TMP_OUTPUT)
R2_FILE = '{0}2.fq'.format(TMP_OUTPUT)


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
    parser.add_argument('--log', default='metaART.log', type=argparse.FileType('w'), help='Log file name')
    parser.add_argument('fasta', metavar='MULTIFASTA',
                        help='Input multi-fasta of all sequences')
    parser.add_argument('output_base', metavar='OUTPUT_BASE',
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

    seq_index = SeqIO.index(args.fasta, 'fasta')
    with open('{0}1.fq'.format(args.output_base), 'w') as output_R1, \
            open('{0}2.fq'.format(args.output_base), 'w') as output_R2:
        try:
            for seq_id in profile:
                coverage = float(profile[seq_id] * args.max_coverage)
                print 'Requesting {0} coverage for {1}'.format(coverage, seq_id)
                ref_seq = seq_index[seq_id]
                ref_len = len(ref_seq)
                SeqIO.write([ref_seq], TMP_INPUT, 'fasta')
                subprocess.call([args.art_path,
                                 '-p',   # paired-end sequencing
                                 '-na',  # no alignment file
                                 '-rs', str(args.seed),
                                 '-m', str(args.insert_len),
                                 '-s', str(args.insert_sd),
                                 '-l', str(args.read_len),
                                 '-f', str(coverage),
                                 '-i', TMP_INPUT,
                                 '-o', TMP_OUTPUT], stdout=args.log, stderr=args.log)

                # count generated reads
                r1_n = 0
                for seq in SeqIO.parse(R1_FILE, 'fastq'):
                    r1_n += 1

                r2_n = 0
                for seq in SeqIO.parse(R1_FILE, 'fastq'):
                    r2_n += 1

                effective_cov = args.read_len * (r1_n + r2_n) / float(ref_len)
                print 'Generated {0} paired-end reads for {1}, {2:.3f} coverage'.format(r1_n, seq_id, effective_cov)
                if r1_n != r2_n:
                    print 'Error: paired-end counts do not match {0} vs {1}'.format(r1_n, r2_n)
                    sys.exit(1)

                with open(R1_FILE, 'r') as tmp_h:
                    output_R1.write(tmp_h.read())
                    os.remove(tmp_h.name)

                with open(R2_FILE, 'r') as tmp_h:
                    output_R2.write(tmp_h.read())
                    os.remove(tmp_h.name)

                os.remove(TMP_INPUT)
        finally:
            if seq_index:
                seq_index.close()
